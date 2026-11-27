#!/usr/bin/env python3
"""
Full GEO Series (GSE) scan → omics_v1.db
- No disease filter; all gse[Entry Type]
- Resumable via sync_cursors (source=geo_full_scan)
- Rate-limited NCBI E-utilities (~3 req/s)
- Does NOT download SOFT/matrix bodies; stores URIs/metadata only
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
import sqlite3
import sys
import time
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Any
from urllib.parse import quote

import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pipeline"))
from paths import data_dir as _data_dir, logs_dir as _logs_dir
_DATA_DIR = _data_dir(); _DATA_DIR.mkdir(parents=True, exist_ok=True)
_LOGS_DIR = _logs_dir(); _LOGS_DIR.mkdir(parents=True, exist_ok=True)
DB_PATH = _DATA_DIR / "omics_v1.db"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "omics-sniffer"
EMAIL = "omics-sniffer@localhost"
TERM = "gse[Entry Type]"
CURSOR_SOURCE = "geo_full_scan"
BATCH_SIZE = 200          # IDs per esearch + esummary cycle
ESUMMARY_CHUNK = 200      # esummary id batch (same as BATCH_SIZE)
MIN_INTERVAL = 0.34       # ~3 req/s
LOG_EVERY = 500           # log every N series processed in this run
MAX_RETRIES = 8
TIMEOUT = 120
TZ = timezone(timedelta(hours=8))  # Asia/Shanghai

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "geo_full_scan.log"),
    ],
)
log = logging.getLogger("geo_full_scan")

_last_req = 0.0


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def throttle() -> None:
    global _last_req
    elapsed = time.monotonic() - _last_req
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_req = time.monotonic()


def ncbi_get(endpoint: str, params: dict) -> Any:
    params = {**params, "tool": TOOL, "email": EMAIL}
    url = f"{BASE}/{endpoint}"
    last_err: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        throttle()
        try:
            r = requests.get(url, params=params, timeout=TIMEOUT)
            if r.status_code == 429 or r.status_code >= 500:
                wait = min(2 ** attempt, 60)
                log.warning("HTTP %s from NCBI (attempt %s), sleep %ss", r.status_code, attempt, wait)
                time.sleep(wait)
                continue
            r.raise_for_status()
            # NCBI sometimes returns error XML/JSON with 200
            ctype = r.headers.get("content-type", "")
            if "json" in ctype or params.get("retmode") == "json":
                data = r.json()
                if isinstance(data, dict) and data.get("error"):
                    raise RuntimeError(f"NCBI JSON error: {data['error']}")
                return data
            return r.text
        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 60)
            log.warning("Request failed attempt %s: %s; sleep %ss", attempt, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"NCBI request failed after {MAX_RETRIES} retries: {last_err}")


def title_norm(title: str | None) -> str:
    if not title:
        return ""
    t = title.lower().strip()
    t = re.sub(r"\s+", " ", t)
    return t


def content_hash(obj: Any) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=False, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:32]


def guess_format(filename: str | None, suppfile: str | None) -> str | None:
    if filename:
        lower = filename.lower()
        for ext in (".csv.gz", ".tsv.gz", ".txt.gz", ".h5ad", ".h5", ".mtx.gz",
                    ".soft.gz", ".xml.gz", ".tar", ".cel.gz", ".fastq.gz",
                    ".csv", ".tsv", ".txt", ".soft", ".cel", ".bw", ".bigwig"):
            if lower.endswith(ext):
                return ext.lstrip(".")
    if suppfile:
        # e.g. "CSV", "CEL", "TXT, CEL"
        return suppfile.split(",")[0].strip().lower() or None
    return None


def artifact_role_from_name(name: str | None, kind_hint: str) -> str:
    n = (name or "").lower()
    if "count" in n:
        return "counts"
    if "tpm" in n:
        return "tpm"
    if "fpkm" in n or "rpkm" in n:
        return "fpkm"
    if "matrix" in n:
        return "expression_matrix"
    if "soft" in n:
        return "soft"
    if "clinical" in n or "phenotype" in n or "meta" in n:
        return "phenotype"
    if kind_hint == "processed":
        return "processed"
    return "other"


def connect() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH, timeout=120)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA cache_size=-200000")  # ~200MB
    return conn


def load_cursor(conn: sqlite3.Connection) -> dict:
    row = conn.execute(
        "SELECT cursor_value, cursor_kind FROM sync_cursors WHERE source=?",
        (CURSOR_SOURCE,),
    ).fetchone()
    if not row or not row[0]:
        return {"retstart": 0, "total": None, "ingested": 0}
    try:
        return json.loads(row[0])
    except json.JSONDecodeError:
        return {"retstart": 0, "total": None, "ingested": 0}


def save_cursor(conn: sqlite3.Connection, cur: dict, error: str | None = None) -> None:
    conn.execute(
        """
        INSERT INTO sync_cursors(source, cursor_value, cursor_kind, last_run_at, last_error)
        VALUES (?, ?, ?, ?, ?)
        ON CONFLICT(source) DO UPDATE SET
            cursor_value=excluded.cursor_value,
            cursor_kind=excluded.cursor_kind,
            last_run_at=excluded.last_run_at,
            last_error=excluded.last_error
        """,
        (CURSOR_SOURCE, json.dumps(cur), "retstart", now_iso(), error),
    )
    conn.commit()


def upsert_series(conn: sqlite3.Connection, doc: dict, fetched_at: str) -> bool:
    """Upsert one GSE into src + unified tables. Returns True if GSE processed."""
    if doc.get("entrytype") != "GSE":
        return False
    if "error" in doc:
        return False

    acc = doc.get("accession") or ""
    if not acc.startswith("GSE"):
        # sometimes accession missing; build from gse field
        gse_num = str(doc.get("gse") or "").strip()
        if gse_num.isdigit():
            acc = f"GSE{gse_num}"
        else:
            return False

    title = doc.get("title") or ""
    summary = doc.get("summary") or ""
    organism = doc.get("taxon") or doc.get("samplestaxa") or None
    gpl = doc.get("gpl") or ""
    platforms = []
    if gpl:
        for p in str(gpl).split(";"):
            p = p.strip()
            if not p:
                continue
            platforms.append(p if p.upper().startswith("GPL") else f"GPL{p}")

    try:
        n_samples = int(doc.get("n_samples") or 0)
    except (TypeError, ValueError):
        n_samples = 0

    samples = doc.get("samples") or []
    if isinstance(samples, dict):
        samples = list(samples.values())
    pubmedids = doc.get("pubmedids") or []
    if isinstance(pubmedids, dict):
        # sometimes {"pubmedid": ["1","2"]} or list-like
        pubmedids = pubmedids.get("pubmedid") or pubmedids.get("PubmedId") or list(pubmedids.values())
    if not isinstance(pubmedids, list):
        pubmedids = [pubmedids] if pubmedids else []
    pubmedids = [str(x).strip() for x in pubmedids if str(x).strip()]

    pdat = doc.get("pdat") or None
    suppfile = doc.get("suppfile") or ""
    ftplink = (doc.get("ftplink") or "").rstrip("/") + ("/" if doc.get("ftplink") else "")
    bioproject = doc.get("bioproject") or ""
    gdstype = doc.get("gdstype") or ""

    # supplementary_files JSON: types + ftp pointers (no file bodies)
    supp_meta = {
        "suppfile_types": [s.strip() for s in str(suppfile).split(",") if s.strip()],
        "ftplink": ftplink or None,
        "bioproject": bioproject or None,
    }
    # Keep a compact raw snapshot (drop huge sample list duplication in raw)
    # raw_xml no longer stored

    hash_payload = {
        "accession": acc,
        "title": title,
        "organism": organism,
        "platforms": platforms,
        "n_samples": n_samples,
        "gdstype": gdstype,
        "pdat": pdat,
        "suppfile": suppfile,
        "pubmedids": pubmedids,
        "ftplink": ftplink,
    }
    ch = content_hash(hash_payload)

    # --- src_geo_series ---
    conn.execute(
        """
        INSERT INTO src_geo_series(
            accession, title, summary, organism, platform, sample_count,
            supplementary_files, raw_xml, fetched_at, updatedate
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(accession) DO UPDATE SET
            title=excluded.title,
            summary=excluded.summary,
            organism=excluded.organism,
            platform=excluded.platform,
            sample_count=excluded.sample_count,
            supplementary_files=excluded.supplementary_files,
            raw_xml=NULL,
            fetched_at=excluded.fetched_at,
            updatedate=excluded.updatedate
        """,
        (
            acc,
            title,
            summary,
            organism,
            ";".join(platforms) if platforms else None,
            n_samples,
            json.dumps(supp_meta, ensure_ascii=False),
            None,  # raw_xml omitted; re-fetch by accession if needed
            fetched_at,
            pdat,
        ),
    )

    canonical_id = f"geo:{acc}"
    # --- datasets ---
    existing = conn.execute(
        "SELECT first_seen_at FROM datasets WHERE canonical_id=?", (canonical_id,)
    ).fetchone()
    first_seen = existing[0] if existing else fetched_at
    access = "public"  # GEO Series listed here are public metadata; controlled noted elsewhere
    conn.execute(
        """
        INSERT INTO datasets(
            canonical_id, title, title_norm, organism, tissue, disease_tags,
            assay, platforms, n_samples, access, license, confidence,
            content_hash, first_seen_at, updated_at, probe_status, last_checked
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(canonical_id) DO UPDATE SET
            title=excluded.title,
            title_norm=excluded.title_norm,
            organism=excluded.organism,
            assay=excluded.assay,
            platforms=excluded.platforms,
            n_samples=excluded.n_samples,
            content_hash=excluded.content_hash,
            updated_at=excluded.updated_at,
            probe_status=excluded.probe_status,
            last_checked=excluded.last_checked
        """,
        (
            canonical_id,
            title,
            title_norm(title),
            organism,
            None,
            None,  # disease_tags not inferred in full scan
            gdstype or None,
            json.dumps(platforms, ensure_ascii=False),
            n_samples,
            access,
            None,
            "metadata_only",
            ch,
            first_seen,
            fetched_at,
            "metadata_fetched",
            fetched_at,
        ),
    )

    # --- dataset_members ---
    geo_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"
    conn.execute(
        """
        INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
        VALUES (?, 'geo', ?, 'series', ?)
        ON CONFLICT(source, accession) DO UPDATE SET
            dataset_id=excluded.dataset_id,
            role=excluded.role,
            external_url=excluded.external_url
        """,
        (canonical_id, acc, geo_url),
    )
    if bioproject:
        # bioproject may be like PRJNA123 or numeric
        bp = str(bioproject).strip()
        if bp and not bp.upper().startswith("PRJ"):
            # leave as-is; NCBI sometimes returns numeric
            pass
        try:
            conn.execute(
                """
                INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
                VALUES (?, 'ncbi', ?, 'bioproject', ?)
                ON CONFLICT(source, accession) DO UPDATE SET
                    dataset_id=excluded.dataset_id,
                    role=excluded.role
                """,
                (
                    canonical_id,
                    bp if bp.upper().startswith("PRJ") else f"bioproject:{bp}",
                    f"https://www.ncbi.nlm.nih.gov/bioproject/{bp}",
                ),
            )
        except sqlite3.IntegrityError:
            pass

    # --- samples (from esummary list; may be truncated vs n_samples) ---
    for s in samples:
        if not isinstance(s, dict):
            continue
        gsm = s.get("accession") or ""
        if not gsm.startswith("GSM"):
            continue
        sample_id = f"geo:{gsm}"
        stitle = s.get("title") or ""
        attrs = {"geo_listed_in_esummary": True, "series": acc}
        conn.execute(
            """
            INSERT INTO samples(
                sample_id, dataset_id, title, organism, tissue, disease,
                condition, sex, age, attributes, run_ids, content_hash, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(sample_id) DO UPDATE SET
                dataset_id=excluded.dataset_id,
                title=excluded.title,
                organism=excluded.organism,
                attributes=excluded.attributes,
                content_hash=excluded.content_hash,
                updated_at=excluded.updated_at
            """,
            (
                sample_id,
                canonical_id,
                stitle,
                organism,
                None,
                None,
                None,
                None,
                None,
                json.dumps(attrs, ensure_ascii=False),
                None,
                content_hash({"gsm": gsm, "title": stitle, "series": acc}),
                fetched_at,
            ),
        )
        conn.execute(
            """
            INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
            VALUES (?, 'geo', ?, 'sample', ?)
            ON CONFLICT(source, accession) DO UPDATE SET
                dataset_id=excluded.dataset_id,
                role=excluded.role
            """,
            (
                canonical_id,
                gsm,
                f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}",
            ),
        )

    # --- artifacts: FTP pointers only ---
    def upsert_artifact(artifact_id: str, kind: str, role: str, title: str,
                        filename: str | None, fmt: str | None, uri: str) -> None:
        conn.execute(
            """
            INSERT INTO dataset_artifacts(
                artifact_id, dataset_id, sample_id, kind, role, title, filename,
                format, uri, size_bytes, access, license, source, derived_from,
                probe_status, last_checked, content_hash
            ) VALUES (?, ?, NULL, ?, ?, ?, ?, ?, ?, NULL, ?, NULL, 'geo', NULL, ?, ?, ?)
            ON CONFLICT(artifact_id) DO UPDATE SET
                uri=excluded.uri,
                format=excluded.format,
                title=excluded.title,
                probe_status=excluded.probe_status,
                last_checked=excluded.last_checked,
                content_hash=excluded.content_hash
            """,
            (
                artifact_id,
                canonical_id,
                kind,
                role,
                title,
                filename,
                fmt,
                uri,
                access,
                "uri_only",
                fetched_at,
                content_hash({"uri": uri, "kind": kind, "role": role}),
            ),
        )

    if ftplink:
        upsert_artifact(
            f"geo:{acc}:ftp_root",
            "raw",
            "ftp_root",
            f"{acc} GEO FTP root",
            None,
            None,
            ftplink,
        )
        # Conventional GEO layout pointers (bodies not downloaded)
        upsert_artifact(
            f"geo:{acc}:soft_dir",
            "supplementary",
            "soft",
            f"{acc} SOFT directory",
            None,
            "soft",
            ftplink + "soft/",
        )
        upsert_artifact(
            f"geo:{acc}:matrix_dir",
            "processed",
            "expression_matrix",
            f"{acc} Series Matrix directory",
            None,
            "matrix",
            ftplink + "matrix/",
        )
        if supp_meta["suppfile_types"]:
            upsert_artifact(
                f"geo:{acc}:suppl_dir",
                "supplementary",
                artifact_role_from_name(suppfile, "supplementary"),
                f"{acc} supplementary files ({suppfile})",
                None,
                guess_format(None, suppfile),
                ftplink + "suppl/",
            )

    # --- publications ---
    for pmid in pubmedids:
        if not pmid.isdigit():
            continue
        conn.execute(
            """
            INSERT INTO publications(pmid, doi, title, authors, journal, year, abstract, url, fetched_at)
            VALUES (?, NULL, NULL, NULL, NULL, NULL, NULL, ?, ?)
            ON CONFLICT(pmid) DO UPDATE SET
                url=COALESCE(publications.url, excluded.url),
                fetched_at=excluded.fetched_at
            """,
            (
                pmid,
                f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                fetched_at,
            ),
        )
        conn.execute(
            """
            INSERT INTO dataset_publications(dataset_id, pmid, relation, source)
            VALUES (?, ?, 'primary', 'geo_link')
            ON CONFLICT(dataset_id, pmid) DO NOTHING
            """,
            (canonical_id, pmid),
        )

    return True


def fetch_total() -> int:
    data = ncbi_get(
        "esearch.fcgi",
        {"db": "gds", "term": TERM, "retmax": 0, "retmode": "json"},
    )
    return int(data["esearchresult"]["count"])


def fetch_id_page(retstart: int, retmax: int) -> list[str]:
    data = ncbi_get(
        "esearch.fcgi",
        {
            "db": "gds",
            "term": TERM,
            "retstart": retstart,
            "retmax": retmax,
            "retmode": "json",
        },
    )
    return list(data["esearchresult"].get("idlist") or [])


def fetch_summaries(ids: list[str]) -> list[dict]:
    if not ids:
        return []
    data = ncbi_get(
        "esummary.fcgi",
        {"db": "gds", "id": ",".join(ids), "retmode": "json"},
    )
    result = data.get("result") or {}
    docs = []
    for uid in result.get("uids") or ids:
        doc = result.get(str(uid))
        if doc and isinstance(doc, dict) and "error" not in doc:
            docs.append(doc)
    return docs


def table_counts(conn: sqlite3.Connection) -> dict[str, int]:
    out = {}
    for t in (
        "datasets",
        "samples",
        "dataset_members",
        "dataset_artifacts",
        "src_geo_series",
        "publications",
        "dataset_publications",
    ):
        out[t] = conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
    return out


def main() -> int:
    log.info("DB=%s", DB_PATH)
    conn = connect()
    cur = load_cursor(conn)
    retstart = int(cur.get("retstart") or 0)
    ingested_prev = int(cur.get("ingested") or 0)

    log.info("Counting total GSE via esearch...")
    total = fetch_total()
    log.info("Total GSE in NCBI GEO: %s", total)
    cur["total"] = total

    if retstart >= total:
        log.info("Cursor retstart=%s already >= total=%s; scan complete.", retstart, total)
        save_cursor(conn, cur, error=None)
        log.info("Counts: %s", table_counts(conn))
        conn.close()
        return 0

    log.info("Resuming at retstart=%s (previously ingested≈%s)", retstart, ingested_prev)
    save_cursor(conn, cur)

    run_ingested = 0
    run_start = time.monotonic()
    errors = 0

    while retstart < total:
        try:
            ids = fetch_id_page(retstart, BATCH_SIZE)
            if not ids:
                log.warning("Empty id page at retstart=%s; advancing by batch", retstart)
                retstart += BATCH_SIZE
                cur["retstart"] = retstart
                save_cursor(conn, cur)
                continue

            docs = fetch_summaries(ids)
            fetched_at = now_iso()
            batch_ok = 0
            conn.execute("BEGIN")
            try:
                for doc in docs:
                    if upsert_series(conn, doc, fetched_at):
                        batch_ok += 1
                conn.commit()
            except Exception:
                conn.rollback()
                raise

            retstart += len(ids)
            run_ingested += batch_ok
            cur["retstart"] = retstart
            cur["ingested"] = ingested_prev + run_ingested
            cur["total"] = total
            save_cursor(conn, cur, error=None)

            if run_ingested % LOG_EVERY < BATCH_SIZE or retstart >= total:
                elapsed = max(time.monotonic() - run_start, 1e-6)
                rate = run_ingested / elapsed
                remaining = max(total - retstart, 0)
                eta_s = remaining / rate if rate > 0 else float("inf")
                counts = table_counts(conn)
                log.info(
                    "progress retstart=%s/%s (+%s this run) rate=%.1f GSE/s ETA≈%.1f min | %s",
                    retstart,
                    total,
                    run_ingested,
                    rate,
                    eta_s / 60.0,
                    counts,
                )
        except Exception as e:
            errors += 1
            log.exception("Batch failed at retstart=%s: %s", retstart, e)
            cur["retstart"] = retstart
            save_cursor(conn, cur, error=str(e)[:500])
            # back off then continue from same retstart
            time.sleep(min(5 * errors, 60))
            if errors >= 50:
                log.error("Too many errors; stopping for resume later.")
                conn.close()
                return 2
            continue

    elapsed = time.monotonic() - run_start
    counts = table_counts(conn)
    log.info(
        "FULL SCAN COMPLETE. retstart=%s total=%s run_ingested=%s elapsed=%.1fs errors=%s",
        retstart,
        total,
        run_ingested,
        elapsed,
        errors,
    )
    log.info("Final counts: %s", counts)
    cur["retstart"] = retstart
    cur["ingested"] = ingested_prev + run_ingested
    cur["total"] = total
    cur["completed_at"] = now_iso()
    save_cursor(conn, cur, error=None)
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
