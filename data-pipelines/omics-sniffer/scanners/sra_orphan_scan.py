#!/usr/bin/env python3
"""
Ingest SRA-only BioProjects (no GEO Series link) into omics_v1.db.

Strategy:
1. Page bioproject db term="Raw sequence reads"[All Fields]
2. Batch elink bioproject→gds; skip if any GEO link
3. Skip if PRJ accession already in dataset_members
4. Create dataset canonical_id=sra:PRJNA… from bioproject esummary
5. Enrich with up to MAX_SRA_IDS_PER_PROJECT experiments via elink+esummary
6. Resumable sync_cursors source=sra_orphan_scan
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

import requests

from http_util import TIMEOUT, http_get, http_post

from paths import data_dir as _data_dir, logs_dir as _logs_dir
_DATA_DIR = _data_dir(); _DATA_DIR.mkdir(parents=True, exist_ok=True)
_LOGS_DIR = _logs_dir(); _LOGS_DIR.mkdir(parents=True, exist_ok=True)
DB_PATH = _DATA_DIR / "omics_v1.db"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "omics-sniffer"
EMAIL = "omics-sniffer@localhost"
CURSOR_SOURCE = "sra_orphan_scan"
TERM = '"Raw sequence reads"[All Fields]'
PAGE_SIZE = 200
ELINK_BATCH = 80
ESUMMARY_BATCH = 200
MAX_SRA_IDS_PER_PROJECT = 100  # cap enrichment per orphan
MIN_INTERVAL = 0.34
MAX_RETRIES = 8
TZ = timezone(timedelta(hours=8))
COMMIT_EVERY_PAGES = 1  # release write lock every page
# If True, still create dataset even when SRA elink fails / empty
CREATE_WITHOUT_SRA = True

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "sra_orphan_scan.log"),
    ],
)
log = logging.getLogger("sra_orphan_scan")
_last_req = 0.0


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def throttle() -> None:
    global _last_req
    elapsed = time.monotonic() - _last_req
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_req = time.monotonic()


def ncbi_request(method: str, endpoint: str, params: dict | None = None, data: Any = None) -> Any:
    params = {**(params or {}), "tool": TOOL, "email": EMAIL}
    url = f"{BASE}/{endpoint}"
    last_err: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        throttle()
        try:
            if method == "GET":
                r = http_get(url, params=params, timeout=TIMEOUT)
            else:
                body = data if data is not None else params
                if isinstance(body, dict):
                    body = {**body, "tool": TOOL, "email": EMAIL}
                    r = http_post(url, data=body, timeout=TIMEOUT)
                else:
                    body = list(body) + [("tool", TOOL), ("email", EMAIL)]
                    r = http_post(url, data=body, timeout=TIMEOUT)
            if r.status_code == 429 or r.status_code >= 500:
                wait = min(2 ** attempt, 60)
                log.warning("HTTP %s (attempt %s), sleep %ss", r.status_code, attempt, wait)
                time.sleep(wait)
                continue
            r.raise_for_status()
            want_json = (params or {}).get("retmode") == "json" or (
                isinstance(data, dict) and data.get("retmode") == "json"
            ) or (
                isinstance(data, list) and any(k == "retmode" and v == "json" for k, v in data)
            )
            ctype = r.headers.get("content-type", "")
            if want_json or "json" in ctype:
                return r.json()
            return r.text
        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 60)
            log.warning("Request failed attempt %s: %s; sleep %ss", attempt, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"NCBI request failed after {MAX_RETRIES}: {last_err}")


def content_hash(obj: Any) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=False, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:32]


def title_norm(title: str | None) -> str:
    if not title:
        return ""
    return re.sub(r"\s+", " ", title.lower().strip())


def connect() -> sqlite3.Connection:
    # timeout= is busy-handler seconds; also set PRAGMA for clarity
    conn = sqlite3.connect(DB_PATH, timeout=120)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA cache_size=-400000")
    conn.execute("PRAGMA busy_timeout=120000")
    return conn


def load_cursor(conn: sqlite3.Connection) -> dict:
    row = conn.execute(
        "SELECT cursor_value FROM sync_cursors WHERE source=?", (CURSOR_SOURCE,)
    ).fetchone()
    if not row or not row[0]:
        return {
            "retstart": 0,
            "total": None,
            "scanned": 0,
            "geo_linked_skipped": 0,
            "already_in_db_skipped": 0,
            "orphans_ingested": 0,
            "experiments_upserted": 0,
            "samples_created": 0,
            "errors": 0,
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
        (CURSOR_SOURCE, json.dumps(cur), "retstart", now_iso(), error),
    )
    conn.commit()


def load_known_bioprojects(conn: sqlite3.Connection) -> set[str]:
    known: set[str] = set()
    for (acc,) in conn.execute(
        "SELECT accession FROM dataset_members WHERE accession LIKE 'PRJ%' OR role='bioproject'"
    ):
        known.add(acc)
        if acc.startswith("bioproject:"):
            known.add(acc.split(":", 1)[1])
    return known


def elink_bioproject_to_gds(uids: list[str]) -> dict[str, list[str]]:
    if not uids:
        return {}
    data = [("dbfrom", "bioproject"), ("db", "gds"), ("retmode", "json")]
    for u in uids:
        data.append(("id", str(u)))
    js = ncbi_request("POST", "elink.fcgi", data=data)
    out: dict[str, list[str]] = {}
    for linkset in js.get("linksets") or []:
        ids = linkset.get("ids") or []
        if not ids:
            continue
        key = str(ids[0])
        gds: list[str] = []
        for ldb in linkset.get("linksetdbs") or []:
            if ldb.get("dbto") == "gds":
                gds.extend(str(x) for x in (ldb.get("links") or []))
        out[key] = gds
    return out


def elink_bioproject_to_sra(uid: str) -> list[str]:
    data = [("dbfrom", "bioproject"), ("db", "sra"), ("retmode", "json"), ("id", str(uid))]
    js = ncbi_request("POST", "elink.fcgi", data=data)
    sra: list[str] = []
    for linkset in js.get("linksets") or []:
        for ldb in linkset.get("linksetdbs") or []:
            if ldb.get("dbto") == "sra":
                sra.extend(str(x) for x in (ldb.get("links") or []))
    return sra


def esummary_bioproject(uids: list[str]) -> dict[str, dict]:
    if not uids:
        return {}
    js = ncbi_request(
        "GET",
        "esummary.fcgi",
        params={"db": "bioproject", "id": ",".join(uids), "retmode": "json"},
    )
    result = js.get("result") or {}
    out = {}
    for uid in uids:
        doc = result.get(str(uid))
        if doc and "error" not in doc:
            out[str(uid)] = doc
    return out


def parse_sra_esummary(uid: str, doc: dict) -> dict | None:
    expxml = doc.get("expxml") or ""
    runs_xml = doc.get("runs") or ""
    m_exp = re.search(r'<Experiment\s+acc="([^"]+)"', expxml)
    if not m_exp:
        return None
    title_m = re.search(r"<Title>(.*?)</Title>", expxml)
    bp_m = re.search(r"<Bioproject>([^<]+)</Bioproject>", expxml)
    bs_m = re.search(r"<Biosample>([^<]+)</Biosample>", expxml)
    study_m = re.search(r'<Study\s+acc="([^"]+)"', expxml)
    srs_m = re.search(r'<Sample\s+acc="([^"]+)"', expxml)
    instrument_m = re.search(r'instrument_model="([^"]+)"', expxml)
    stats_m = re.search(
        r'<Statistics[^>]*total_spots="(\d+)"[^>]*total_bases="(\d+)"', expxml
    )
    lib_m = re.search(r"<LIBRARY_STRATEGY>([^<]+)</LIBRARY_STRATEGY>", expxml)
    org_m = re.search(r'ScientificName="([^"]+)"', expxml)
    run_accs = re.findall(r'<Run\s+acc="([^"]+)"', runs_xml)
    return {
        "uid": uid,
        "accession": m_exp.group(1),
        "title": title_m.group(1) if title_m else m_exp.group(1),
        "bioproject": bp_m.group(1).strip() if bp_m else None,
        "biosample": bs_m.group(1).strip() if bs_m else None,
        "srs": srs_m.group(1) if srs_m else None,
        "study": study_m.group(1) if study_m else None,
        "platform": instrument_m.group(1) if instrument_m else None,
        "spots": int(stats_m.group(1)) if stats_m else None,
        "bases": int(stats_m.group(2)) if stats_m else None,
        "run_accessions": run_accs,
        "library_strategy": lib_m.group(1) if lib_m else None,
        "organism": org_m.group(1) if org_m else None,
    }


def fetch_sra_parsed(uids: list[str]) -> list[dict]:
    out: list[dict] = []
    for i in range(0, len(uids), ESUMMARY_BATCH):
        chunk = uids[i : i + ESUMMARY_BATCH]
        js = ncbi_request(
            "GET",
            "esummary.fcgi",
            params={"db": "sra", "id": ",".join(chunk), "retmode": "json"},
        )
        result = js.get("result") or {}
        for uid in chunk:
            doc = result.get(str(uid))
            if not doc or "error" in doc:
                continue
            p = parse_sra_esummary(str(uid), doc)
            if p:
                out.append(p)
    return out


def bp_accession(doc: dict, uid: str) -> str:
    acc = (doc.get("project_acc") or doc.get("accession") or "").strip()
    if acc:
        return acc
    return f"bioproject_uid:{uid}"


def ingest_orphan(conn: sqlite3.Connection, bp_uid: str, bp_doc: dict, known: set[str], cur: dict) -> bool:
    acc = bp_accession(bp_doc, bp_uid)
    if acc in known:
        cur["already_in_db_skipped"] = cur.get("already_in_db_skipped", 0) + 1
        return False

    title = (bp_doc.get("project_title") or bp_doc.get("title") or acc).strip()
    organism = bp_doc.get("organism_name") or bp_doc.get("taxname") or None
    if isinstance(organism, list):
        organism = "; ".join(str(x) for x in organism)
    fetched_at = now_iso()
    canonical_id = f"sra:{acc}" if not acc.startswith("bioproject_uid:") else f"sra:uid:{bp_uid}"

    # Never hold a write transaction across network I/O
    try:
        conn.commit()
    except Exception:
        pass

    # Enrich with capped SRA experiments
    sra_ids: list[str] = []
    truncated = False
    try:
        sra_ids = elink_bioproject_to_sra(bp_uid)
    except Exception as e:
        log.warning("elink sra failed for %s: %s", acc, e)
    if len(sra_ids) > MAX_SRA_IDS_PER_PROJECT:
        sra_ids = sra_ids[:MAX_SRA_IDS_PER_PROJECT]
        truncated = True

    parsed: list[dict] = []
    if sra_ids:
        try:
            parsed = fetch_sra_parsed(sra_ids)
        except Exception as e:
            cur["errors"] = cur.get("errors", 0) + 1
            log.error("sra esummary failed for %s: %s", acc, e)

    if not parsed and not CREATE_WITHOUT_SRA:
        return False

    strategies = sorted({p.get("library_strategy") or "" for p in parsed if p.get("library_strategy")})
    platforms = sorted({p.get("platform") or "" for p in parsed if p.get("platform")})
    orgs = sorted({p.get("organism") or "" for p in parsed if p.get("organism")})
    if orgs and not organism:
        organism = orgs[0] if len(orgs) == 1 else "; ".join(orgs[:5])

    # samples by biosample
    sample_map: dict[str, dict] = {}
    for p in parsed:
        key = p.get("biosample") or p.get("srs") or p["accession"]
        if key not in sample_map:
            sample_map[key] = {
                "biosample": p.get("biosample"),
                "srs": p.get("srs"),
                "title": p.get("title"),
                "organism": p.get("organism"),
                "runs": set(),
            }
        sample_map[key]["runs"].update(p.get("run_accessions") or [])

    n_samples = len(sample_map) or None
    ch = content_hash({"acc": acc, "title": title, "n": n_samples, "trunc": truncated})
    existing = conn.execute(
        "SELECT first_seen_at FROM datasets WHERE canonical_id=?", (canonical_id,)
    ).fetchone()
    first_seen = existing[0] if existing else fetched_at

    conn.execute(
        """
        INSERT INTO datasets(
            canonical_id, title, title_norm, organism, tissue, disease_tags,
            assay, platforms, n_samples, access, license, confidence,
            content_hash, first_seen_at, updated_at, probe_status, last_checked
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(canonical_id) DO UPDATE SET
            title=excluded.title, title_norm=excluded.title_norm, organism=excluded.organism,
            assay=excluded.assay, platforms=excluded.platforms, n_samples=excluded.n_samples,
            content_hash=excluded.content_hash, updated_at=excluded.updated_at,
            probe_status=excluded.probe_status, last_checked=excluded.last_checked
        """,
        (
            canonical_id, title, title_norm(title), organism, None, None,
            "; ".join(strategies) if strategies else None,
            json.dumps(platforms),
            n_samples, "public", None, "metadata_only", ch, first_seen, fetched_at,
            "metadata_fetched" + ("_truncated" if truncated else ""), fetched_at,
        ),
    )

    if acc.startswith("PRJ"):
        try:
            conn.execute(
                """
                INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
                VALUES (?, 'sra', ?, 'bioproject', ?)
                ON CONFLICT(source, accession) DO UPDATE SET
                    dataset_id=excluded.dataset_id, role=excluded.role
                """,
                (canonical_id, acc, f"https://www.ncbi.nlm.nih.gov/bioproject/{acc}"),
            )
        except sqlite3.IntegrityError:
            pass
        known.add(acc)

    studies = sorted({p.get("study") for p in parsed if p.get("study")})
    for study in studies:
        try:
            conn.execute(
                """
                INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
                VALUES (?, 'sra', ?, 'study', ?)
                ON CONFLICT(source, accession) DO UPDATE SET
                    dataset_id=excluded.dataset_id, role=excluded.role
                """,
                (canonical_id, study, f"https://www.ncbi.nlm.nih.gov/sra?term={study}"),
            )
        except sqlite3.IntegrityError:
            pass

    for p in parsed:
        exp = p["accession"]
        conn.execute(
            """
            INSERT INTO src_sra_experiment(
                accession, bioproject, biosample, title, platform, spots, bases,
                run_accessions, ftp_paths, raw_xml, fetched_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(accession) DO UPDATE SET
                bioproject=excluded.bioproject, biosample=excluded.biosample,
                title=excluded.title, platform=excluded.platform, spots=excluded.spots,
                bases=excluded.bases, run_accessions=excluded.run_accessions,
                raw_xml=NULL, fetched_at=excluded.fetched_at
            """,
            (
                exp, p.get("bioproject") or acc, p.get("biosample"), p.get("title"),
                p.get("platform"), p.get("spots"), p.get("bases"),
                json.dumps(p.get("run_accessions") or []), None,
                None,  # raw_xml omitted
                fetched_at,
            ),
        )
        cur["experiments_upserted"] = cur.get("experiments_upserted", 0) + 1
        try:
            conn.execute(
                """
                INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
                VALUES (?, 'sra', ?, 'experiment', ?)
                ON CONFLICT(source, accession) DO UPDATE SET
                    dataset_id=excluded.dataset_id, role=excluded.role
                """,
                (canonical_id, exp, f"https://www.ncbi.nlm.nih.gov/sra/{exp}"),
            )
        except sqlite3.IntegrityError:
            pass

    for key, meta in sample_map.items():
        bsa = meta["biosample"]
        if bsa:
            sample_id = f"sra:{bsa}"
        elif meta["srs"]:
            sample_id = f"sra:{meta['srs']}"
        else:
            safe = re.sub(r"[^A-Za-z0-9._-]+", "_", key)[:80]
            sample_id = f"sra:{acc}:{safe}"
        attrs = {"biosample": bsa, "srs": meta["srs"], "bioproject": acc, "orphan_sra": True}
        runs = sorted(meta["runs"])
        conn.execute(
            """
            INSERT INTO samples(
                sample_id, dataset_id, title, organism, tissue, disease,
                condition, sex, age, attributes, run_ids, content_hash, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(sample_id) DO UPDATE SET
                dataset_id=excluded.dataset_id, title=excluded.title,
                organism=excluded.organism, attributes=excluded.attributes,
                run_ids=excluded.run_ids, content_hash=excluded.content_hash,
                updated_at=excluded.updated_at
            """,
            (
                sample_id, canonical_id, meta["title"] or key, meta["organism"] or organism,
                None, None, None, None, None, json.dumps(attrs), json.dumps(runs),
                content_hash(attrs), fetched_at,
            ),
        )
        cur["samples_created"] = cur.get("samples_created", 0) + 1

    cur["orphans_ingested"] = cur.get("orphans_ingested", 0) + 1
    known.add(acc)
    conn.commit()  # short write txn; release lock for AE merge
    return True


def main() -> int:
    conn = connect()
    cur = load_cursor(conn)
    known = load_known_bioprojects(conn)
    log.info("Known bioproject accessions in DB: %s", len(known))

    # Re-count if missing or bogus completed-with-total=0 from bad TERM
    if cur.get("total") is None or int(cur.get("total") or 0) == 0:
        if cur.get("completed_at") and int(cur.get("total") or 0) == 0:
            log.warning("Resetting bogus completed cursor (total=0)")
            cur.pop("completed_at", None)
            cur["retstart"] = 0
        js = ncbi_request(
            "GET", "esearch.fcgi",
            params={"db": "bioproject", "term": TERM, "retmode": "json", "retmax": 0},
        )
        cur["total"] = int(js["esearchresult"]["count"])
        save_cursor(conn, cur)
        log.info("Discovery TERM=%s total=%s", TERM, cur["total"])

    total = int(cur["total"])
    retstart = int(cur.get("retstart") or 0)
    if cur.get("completed_at") and retstart >= total and total > 0:
        log.info("Already completed at %s; nothing to do. %s", cur["completed_at"], json.dumps(cur))
        return 0
    log.info("Start sra_orphan_scan TERM=%s retstart=%s total=%s", TERM, retstart, total)

    pages_since_commit = 0
    try:
        while retstart < total:
            try:
                conn.commit()  # flush before page-level NCBI calls
            except Exception:
                pass
            js = ncbi_request(
                "GET", "esearch.fcgi",
                params={
                    "db": "bioproject", "term": TERM, "retmode": "json",
                    "retstart": retstart, "retmax": PAGE_SIZE,
                },
            )
            idlist = js["esearchresult"].get("idlist") or []
            if not idlist:
                log.warning("Empty idlist at retstart=%s; stop", retstart)
                break

            try:
                docs = esummary_bioproject(idlist)
            except Exception as e:
                cur["errors"] = cur.get("errors", 0) + 1
                log.error("esummary failed at %s: %s", retstart, e)
                retstart += len(idlist)
                cur["retstart"] = retstart
                save_cursor(conn, cur, error=str(e))
                continue

            geo_linked: set[str] = set()
            for j in range(0, len(idlist), ELINK_BATCH):
                sub = idlist[j : j + ELINK_BATCH]
                try:
                    link_map = elink_bioproject_to_gds(sub)
                except Exception as e:
                    cur["errors"] = cur.get("errors", 0) + 1
                    log.error("elink gds failed: %s", e)
                    link_map = {}
                for uid, gds_ids in link_map.items():
                    if gds_ids:
                        geo_linked.add(uid)

            for uid in idlist:
                cur["scanned"] = cur.get("scanned", 0) + 1
                doc = docs.get(uid) or {}
                acc = bp_accession(doc, uid) if doc else ""
                if uid in geo_linked:
                    cur["geo_linked_skipped"] = cur.get("geo_linked_skipped", 0) + 1
                    if acc:
                        known.add(acc)
                    continue
                if not doc:
                    continue
                if acc and acc in known:
                    cur["already_in_db_skipped"] = cur.get("already_in_db_skipped", 0) + 1
                    continue
                try:
                    ingest_orphan(conn, uid, doc, known, cur)
                except Exception as e:
                    cur["errors"] = cur.get("errors", 0) + 1
                    log.error("ingest orphan uid=%s failed: %s", uid, e)

            retstart += len(idlist)
            cur["retstart"] = retstart
            pages_since_commit += 1
            if pages_since_commit >= COMMIT_EVERY_PAGES:
                save_cursor(conn, cur)
                pages_since_commit = 0
                log.info(
                    "Progress retstart=%s/%s scanned=%s orphans=%s geo_skip=%s db_skip=%s exp=%s samples=%s errors=%s",
                    retstart, total, cur.get("scanned"), cur.get("orphans_ingested"),
                    cur.get("geo_linked_skipped"), cur.get("already_in_db_skipped"),
                    cur.get("experiments_upserted"), cur.get("samples_created"),
                    cur.get("errors"),
                )

        cur["retstart"] = total
        cur["completed_at"] = now_iso()
        save_cursor(conn, cur)
        log.info("COMPLETED sra_orphan_scan %s", json.dumps(cur))
        return 0
    except KeyboardInterrupt:
        save_cursor(conn, cur, error="interrupted")
        return 130
    except Exception as e:
        save_cursor(conn, cur, error=str(e))
        log.exception("Fatal: %s", e)
        return 1
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
