#!/usr/bin/env python3
"""
MassIVE (UCSD / ProteomeXchange partner) full ingest (scan / staging phase).

Writes ONLY to omics_massive_staging.db — never to omics_v1.db.
Uses datasets_json.jsp catalog (cached locally) for lean metadata + FTP pointers.
Resumable via sync_cursors (source=massive_full_scan).
"""

from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)

import hashlib
import json
import logging
import re
import sqlite3
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from http_util import TIMEOUT, http_get

from paths import CODE_ROOT as ROOT, SCHEMA_PATH, data_dir as _data_dir, logs_dir as _logs_dir
_DATA_DIR = _data_dir(); _DATA_DIR.mkdir(parents=True, exist_ok=True)
_LOGS_DIR = _logs_dir(); _LOGS_DIR.mkdir(parents=True, exist_ok=True)
STAGING_DB = _DATA_DIR / "omics_massive_staging.db"
SCHEMA = SCHEMA_PATH
CATALOG_CACHE = ROOT / "massive_catalog_cache.json"
CATALOG_URL = "https://massive.ucsd.edu/ProteoSAFe/datasets_json.jsp"
CURSOR_SOURCE = "massive_full_scan"
MIN_INTERVAL = 0.25
MAX_RETRIES = 8
TZ = timezone(timedelta(hours=8))
SAVE_EVERY = 100
LOG_EVERY = 200
BATCH_COMMIT = 50
DESC_MAX = 4000
CACHE_MAX_AGE_H = 48

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "massive_full_scan.log"),
    ],
)
log = logging.getLogger("massive_full_scan")
_last_req = 0.0

MSV_RE = re.compile(r"^(MSV\d+)$", re.I)
TAXON_RE = re.compile(r"^(.*?)\s*\(NCBITaxon:(\d+)\)\s*$")


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def throttle() -> None:
    global _last_req
    elapsed = time.monotonic() - _last_req
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_req = time.monotonic()


def massive_get(url: str, timeout: Any = TIMEOUT) -> Any:
    last_err: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        throttle()
        try:
            r = http_get(url, timeout=timeout)
            if r.status_code == 429 or r.status_code >= 500:
                wait = min(2**attempt, 30)
                log.warning("HTTP %s %s (attempt %s), sleep %ss", r.status_code, url, attempt, wait)
                time.sleep(wait)
                continue
            if r.status_code in (400, 404):
                return None
            r.raise_for_status()
            return r
        except Exception as e:
            last_err = e
            wait = min(2**attempt, 30)
            log.warning("Request failed attempt %s %s: %s; sleep %ss", attempt, url, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"MassIVE request failed after {MAX_RETRIES}: {last_err}")


def content_hash(obj: Any) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=False, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:32]


def title_norm(title: str | None) -> str:
    if not title:
        return ""
    return re.sub(r"\s+", " ", title.lower().strip())


def connect(path: Path, busy_ms: int = 15000) -> sqlite3.Connection:
    conn = sqlite3.connect(path, timeout=max(busy_ms / 1000, 5))
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute(f"PRAGMA busy_timeout={busy_ms}")
    return conn


def init_staging() -> None:
    need = not STAGING_DB.exists()
    conn = connect(STAGING_DB, busy_ms=5000)
    if need:
        conn.executescript(SCHEMA.read_text())
        conn.commit()
        log.info("Created staging DB %s from schema_v1.sql", STAGING_DB)
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS src_massive_dataset (
            accession              TEXT PRIMARY KEY,  -- MSV######
            title                  TEXT,
            description            TEXT,
            species                TEXT,
            instrument             TEXT,
            keywords               TEXT,
            pi                     TEXT,
            user_name              TEXT,
            px_accession           TEXT,
            file_count             INTEGER,
            file_size_kb           INTEGER,
            status                 TEXT,
            flowname               TEXT,
            created                TEXT,
            created_millis         TEXT,
            quant_analysis         TEXT,
            spectra                TEXT,
            peptides               TEXT,
            proteins               TEXT,
            task                   TEXT,
            site                   TEXT,
            raw_json               TEXT,
            fetched_at             TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_src_massive_fetched_at ON src_massive_dataset(fetched_at);
        CREATE INDEX IF NOT EXISTS idx_src_massive_px ON src_massive_dataset(px_accession);
        CREATE TABLE IF NOT EXISTS massive_merge_state (
            accession TEXT PRIMARY KEY,
            merged_at TEXT
        );
        """
    )
    conn.commit()
    conn.close()


def load_cursor(conn: sqlite3.Connection) -> dict:
    row = conn.execute(
        "SELECT cursor_value FROM sync_cursors WHERE source=?", (CURSOR_SOURCE,)
    ).fetchone()
    if not row or not row[0]:
        return {
            "index": 0,
            "total": None,
            "scanned": 0,
            "skipped_existing": 0,
            "skipped_private": 0,
            "ingested_new": 0,
            "src_upserted": 0,
            "artifacts_created": 0,
            "errors": 0,
            "catalog_fetched_at": None,
        }
    return json.loads(row[0])


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
        (CURSOR_SOURCE, json.dumps(cur), "index", now_iso(), error),
    )
    conn.commit()


def as_int(val: Any) -> int | None:
    if val is None or val == "":
        return None
    try:
        return int(float(str(val).replace(",", "")))
    except (TypeError, ValueError):
        return None



def compact_pi(val: Any) -> str | None:
    """Normalize PI field; drop emails from stored snapshot."""
    if val is None or val == "":
        return None
    if isinstance(val, str):
        return val
    names: list[str] = []
    if isinstance(val, list):
        for item in val:
            if isinstance(item, dict):
                name = item.get("name") or item.get("pi")
                if name:
                    names.append(str(name))
            elif item:
                names.append(str(item))
    elif isinstance(val, dict):
        name = val.get("name") or val.get("pi")
        if name:
            names.append(str(name))
    return "; ".join(dict.fromkeys(names)) if names else None


def parse_species(raw: str | None) -> str | None:
    if not raw:
        return None
    parts = []
    for chunk in str(raw).split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        m = TAXON_RE.match(chunk)
        parts.append(m.group(1).strip() if m else chunk)
    return "; ".join(dict.fromkeys(parts)) if parts else None


def load_catalog(force: bool = False) -> list[dict]:
    use_cache = False
    if CATALOG_CACHE.exists() and not force:
        age_h = (time.time() - CATALOG_CACHE.stat().st_mtime) / 3600
        if age_h < CACHE_MAX_AGE_H:
            use_cache = True
    if use_cache:
        log.info("Loading MassIVE catalog cache %s", CATALOG_CACHE)
        data = json.loads(CATALOG_CACHE.read_text())
        ds = data.get("datasets") if isinstance(data, dict) else data
        if not isinstance(ds, list):
            raise RuntimeError("Bad MassIVE catalog cache")
        log.info("Cache has %s datasets (age ok)", len(ds))
        return ds

    log.info("Fetching MassIVE catalog %s ...", CATALOG_URL)
    r = massive_get(CATALOG_URL, timeout=(15, 180))
    if r is None:
        raise RuntimeError("MassIVE catalog returned empty")
    text = r.text
    start = text.find("{")
    if start < 0:
        raise RuntimeError("MassIVE catalog JSON object not found")
    payload = json.loads(text[start:])
    ds = payload.get("datasets") or []
    if not isinstance(ds, list) or not ds:
        raise RuntimeError("MassIVE catalog datasets empty")
    # lean cache: drop unused heavy-ish keys we don't need later
    lean = []
    keep = {
        "dataset", "title", "description", "species", "instrument", "keywords",
        "pi", "user", "px", "fileCount", "fileSizeKB", "status", "flowname",
        "created", "createdMillis", "quant_analysis", "spectra", "peptides",
        "proteins", "task", "site", "private", "complete", "hash", "id",
        "datasetNum", "modification", "variants", "psms",
    }
    for d in ds:
        if not isinstance(d, dict):
            continue
        lean.append({k: d.get(k) for k in keep if k in d})
    CATALOG_CACHE.write_text(json.dumps({"datasets": lean}, ensure_ascii=False))
    log.info("Cached %s datasets -> %s (%.1f MB)", len(lean), CATALOG_CACHE, CATALOG_CACHE.stat().st_size / 1e6)
    return lean


def staging_ingest(conn: sqlite3.Connection, row: dict, cur: dict) -> None:
    accession = (row.get("dataset") or "").strip().upper()
    if not MSV_RE.match(accession):
        cur["errors"] = cur.get("errors", 0) + 1
        log.warning("skip non-MSV accession: %r", row.get("dataset"))
        return
    if str(row.get("private") or "").lower() in ("1", "true", "yes"):
        cur["skipped_private"] = cur.get("skipped_private", 0) + 1
        return

    exists = conn.execute(
        "SELECT 1 FROM src_massive_dataset WHERE accession=? LIMIT 1", (accession,)
    ).fetchone()
    if exists:
        cur["skipped_existing"] = cur.get("skipped_existing", 0) + 1
        return

    fetched_at = now_iso()
    title = row.get("title") or accession
    desc = row.get("description") or ""
    if isinstance(desc, str) and len(desc) > DESC_MAX:
        desc = desc[:DESC_MAX] + "…"
    species = parse_species(row.get("species"))
    instrument = row.get("instrument")
    if isinstance(instrument, list):
        instrument = "; ".join(str(x) for x in instrument if x)
    keywords = row.get("keywords")
    if isinstance(keywords, list):
        keywords = "; ".join(str(x) for x in keywords if x)
    elif isinstance(keywords, dict):
        keywords = json.dumps(keywords, ensure_ascii=False)
    pi_str = compact_pi(row.get("pi"))
    px = (row.get("px") or "").strip().upper() or None
    if px and not px.startswith("PXD"):
        px = None
    file_count = as_int(row.get("fileCount"))
    file_size_kb = as_int(row.get("fileSizeKB"))

    raw = {
        k: row.get(k)
        for k in (
            "dataset", "datasetNum", "title", "species", "instrument", "keywords",
            "pi", "user", "px", "fileCount", "fileSizeKB", "status", "flowname",
            "created", "createdMillis", "quant_analysis", "spectra", "peptides",
            "proteins", "task", "site", "complete", "hash", "id",
        )
    }
    raw["description"] = desc
    raw["pi"] = pi_str

    conn.execute(
        """
        INSERT INTO src_massive_dataset(
            accession, title, description, species, instrument, keywords, pi,
            user_name, px_accession, file_count, file_size_kb, status, flowname,
            created, created_millis, quant_analysis, spectra, peptides, proteins,
            task, site, raw_json, fetched_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(accession) DO UPDATE SET
            title=excluded.title, description=excluded.description,
            species=excluded.species, instrument=excluded.instrument,
            keywords=excluded.keywords, pi=excluded.pi, user_name=excluded.user_name,
            px_accession=excluded.px_accession, file_count=excluded.file_count,
            file_size_kb=excluded.file_size_kb, status=excluded.status,
            flowname=excluded.flowname, created=excluded.created,
            created_millis=excluded.created_millis,
            quant_analysis=excluded.quant_analysis, spectra=excluded.spectra,
            peptides=excluded.peptides, proteins=excluded.proteins,
            task=excluded.task, site=excluded.site, raw_json=excluded.raw_json,
            fetched_at=excluded.fetched_at
        """,
        (
            accession,
            title,
            desc,
            species,
            instrument,
            keywords,
            pi_str,
            row.get("user"),
            px,
            file_count,
            file_size_kb,
            row.get("status"),
            row.get("flowname"),
            row.get("created"),
            str(row.get("createdMillis") or "") or None,
            row.get("quant_analysis"),
            row.get("spectra"),
            row.get("peptides"),
            row.get("proteins"),
            row.get("task"),
            row.get("site"),
            json.dumps(raw, ensure_ascii=False, default=str),
            fetched_at,
        ),
    )
    cur["src_upserted"] = cur.get("src_upserted", 0) + 1

    canonical_id = f"massive:{accession}"
    ch = content_hash(
        {"acc": accession, "title": title, "species": species, "nf": file_count, "px": px}
    )
    existing = conn.execute(
        "SELECT first_seen_at FROM datasets WHERE canonical_id=?", (canonical_id,)
    ).fetchone()
    first_seen = existing[0] if existing else fetched_at
    platforms = json.dumps([instrument], ensure_ascii=False) if instrument else None

    conn.execute(
        """
        INSERT INTO datasets(
            canonical_id, title, title_norm, organism, tissue, disease_tags,
            assay, platforms, n_samples, access, license, confidence,
            content_hash, first_seen_at, updated_at, probe_status, last_checked
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(canonical_id) DO UPDATE SET
            title=excluded.title, title_norm=excluded.title_norm,
            organism=COALESCE(excluded.organism, datasets.organism),
            assay=COALESCE(excluded.assay, datasets.assay),
            platforms=COALESCE(excluded.platforms, datasets.platforms),
            content_hash=excluded.content_hash, updated_at=excluded.updated_at,
            probe_status=excluded.probe_status, last_checked=excluded.last_checked
        """,
        (
            canonical_id,
            title,
            title_norm(title),
            species,
            None,
            None,
            "proteomics",
            platforms,
            None,
            "public",
            None,
            "metadata_only",
            ch,
            first_seen,
            fetched_at,
            "staging",
            fetched_at,
        ),
    )
    if not existing:
        cur["ingested_new"] = cur.get("ingested_new", 0) + 1

    page_url = f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task={row.get('task')}" if row.get("task") else (
        f"https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id={accession}"
    )
    conn.execute(
        """
        INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
        VALUES (?, 'massive', ?, 'dataset', ?)
        ON CONFLICT(source, accession) DO UPDATE SET
            dataset_id=excluded.dataset_id, role=excluded.role, external_url=excluded.external_url
        """,
        (canonical_id, accession, page_url),
    )
    if px:
        conn.execute(
            """
            INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
            VALUES (?, 'proteomexchange', ?, 'px_link', ?)
            ON CONFLICT(source, accession) DO NOTHING
            """,
            (
                canonical_id,
                px,
                f"http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={px}",
            ),
        )

    conn.execute("DELETE FROM dataset_artifacts WHERE dataset_id=?", (canonical_id,))
    ftp_uri = f"ftp://massive.ucsd.edu/v01/{accession}"
    size_bytes = file_size_kb * 1024 if file_size_kb is not None else None
    for kind, role, uri, title_a, fmt in (
        ("raw", "dataset_ftp", ftp_uri, f"{accession} FTP root", "ftp"),
        ("supplementary", "dataset_page", page_url, f"{accession} MassIVE page", "html"),
    ):
        artifact_id = f"massive:{accession}:{hashlib.sha1(uri.encode()).hexdigest()[:16]}"
        ach = content_hash({"uri": uri, "nf": file_count, "sz": size_bytes})
        conn.execute(
            """
            INSERT INTO dataset_artifacts(
                artifact_id, dataset_id, sample_id, kind, role, title, filename,
                format, uri, size_bytes, access, license, source, derived_from,
                probe_status, last_checked, content_hash
            ) VALUES (?, ?, NULL, ?, ?, ?, ?, ?, ?, ?, 'public', NULL, 'massive', NULL,
                      'listed', ?, ?)
            ON CONFLICT(artifact_id) DO UPDATE SET
                uri=excluded.uri, size_bytes=excluded.size_bytes,
                last_checked=excluded.last_checked, content_hash=excluded.content_hash
            """,
            (
                artifact_id,
                canonical_id,
                kind,
                role,
                title_a,
                accession if kind == "raw" else None,
                fmt,
                uri,
                size_bytes if kind == "raw" else None,
                fetched_at,
                ach,
            ),
        )
        cur["artifacts_created"] = cur.get("artifacts_created", 0) + 1


def main() -> int:
    log.info("Starting MassIVE SCAN-ONLY (staging=%s)", STAGING_DB)
    init_staging()
    st = connect(STAGING_DB, busy_ms=5000)
    cur = load_cursor(st)

    try:
        catalog = load_catalog(force=False)
        total = len(catalog)
        cur["total"] = total
        cur["catalog_fetched_at"] = now_iso()
        idx = int(cur.get("index") or 0)
        if idx < 0:
            idx = 0
        if idx > total:
            idx = total
        save_cursor(st, cur)
        log.info("Resume index=%s / total=%s", idx, total)

        since_log = 0
        since_save = 0
        pending_commit = 0

        while idx < total:
            row = catalog[idx]
            try:
                staging_ingest(st, row if isinstance(row, dict) else {}, cur)
                pending_commit += 1
            except Exception as e:
                cur["errors"] = cur.get("errors", 0) + 1
                acc = (row or {}).get("dataset") if isinstance(row, dict) else None
                log.exception("Failed staging %s: %s", acc, e)
                save_cursor(st, cur, error=str(e)[:500])
                pending_commit = 0

            cur["scanned"] = cur.get("scanned", 0) + 1
            idx += 1
            cur["index"] = idx
            since_log += 1
            since_save += 1

            if pending_commit >= BATCH_COMMIT:
                st.commit()
                pending_commit = 0

            if since_save >= SAVE_EVERY:
                if pending_commit:
                    st.commit()
                    pending_commit = 0
                save_cursor(st, cur)
                since_save = 0

            if since_log >= LOG_EVERY:
                since_log = 0
                log.info(
                    "progress index=%s/%s scanned=%s new=%s skipped=%s src=%s artifacts=%s errors=%s",
                    idx,
                    total,
                    cur.get("scanned"),
                    cur.get("ingested_new"),
                    cur.get("skipped_existing"),
                    cur.get("src_upserted"),
                    cur.get("artifacts_created"),
                    cur.get("errors"),
                )

        if pending_commit:
            st.commit()
        cur["index"] = total
        cur["completed_at"] = now_iso()
        save_cursor(st, cur)
        log.info("Finished MassIVE scan: %s", json.dumps(cur))
        log.info("Staging only — no merge to omics_v1.db")
    except KeyboardInterrupt:
        log.warning("Interrupted at index=%s", cur.get("index"))
        save_cursor(st, cur, error="interrupted")
        st.close()
        return 130
    except Exception as e:
        log.exception("Fatal: %s", e)
        save_cursor(st, cur, error=str(e)[:500])
        st.close()
        return 1

    st.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
