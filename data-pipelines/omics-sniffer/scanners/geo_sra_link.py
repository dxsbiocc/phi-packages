#!/usr/bin/env python3
"""
GEO → SRA association fill for omics_v1.db

- Walks all GSE in src_geo_series
- GDS UID = 200000000 + GSE number; batch elink gds→sra (by-id)
- Batch esummary for linked SRA experiment UIDs (no runinfo / no FASTQ)
- Upserts src_sra_experiment, dataset_members (source=sra), samples.run_ids
- Resumable via sync_cursors source=geo_sra_link
- Rate limit ~3 req/s
"""

from __future__ import annotations

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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pipeline"))
from paths import data_dir as _data_dir, logs_dir as _logs_dir
_DATA_DIR = _data_dir(); _DATA_DIR.mkdir(parents=True, exist_ok=True)
_LOGS_DIR = _logs_dir(); _LOGS_DIR.mkdir(parents=True, exist_ok=True)
DB_PATH = _DATA_DIR / "omics_v1.db"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "omics-sniffer"
EMAIL = "omics-sniffer@localhost"
CURSOR_SOURCE = "geo_sra_link"
ELINK_BATCH = 100
ESUMMARY_BATCH = 200
MAX_SRA_PER_GSE = 500  # cap per series to avoid multi-hour single GSE
MIN_INTERVAL = 0.34
MAX_RETRIES = 8
TIMEOUT = 180
TZ = timezone(timedelta(hours=8))
COMMIT_EVERY = 5
LOG_EVERY = 500

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "geo_sra_link.log"),
    ],
)
log = logging.getLogger("geo_sra_link")
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
                r = requests.get(url, params=params, timeout=TIMEOUT)
            else:
                body = data if data is not None else params
                if isinstance(body, dict):
                    body = {**body, "tool": TOOL, "email": EMAIL}
                    r = requests.post(url, data=body, timeout=TIMEOUT)
                else:
                    body = list(body) + [("tool", TOOL), ("email", EMAIL)]
                    r = requests.post(url, data=body, timeout=TIMEOUT)
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


def connect() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH, timeout=300)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA cache_size=-400000")
    return conn


def load_cursor(conn: sqlite3.Connection) -> dict:
    row = conn.execute(
        "SELECT cursor_value FROM sync_cursors WHERE source=?", (CURSOR_SOURCE,)
    ).fetchone()
    if not row or not row[0]:
        return {
            "offset": 0,
            "gse_total": None,
            "gse_processed": 0,
            "gse_with_sra": 0,
            "experiments_upserted": 0,
            "members_added": 0,
            "samples_updated": 0,
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
        (CURSOR_SOURCE, json.dumps(cur), "gse_offset", now_iso(), error),
    )
    conn.commit()


def gds_uid_for_gse(gse: str) -> int | None:
    m = re.fullmatch(r"GSE(\d+)", gse)
    if not m:
        return None
    return 200_000_000 + int(m.group(1))


def elink_gds_to_sra(gds_uids: list[int]) -> dict[str, list[str]]:
    if not gds_uids:
        return {}
    data = [("dbfrom", "gds"), ("db", "sra"), ("retmode", "json")]
    for u in gds_uids:
        data.append(("id", str(u)))
    js = ncbi_request("POST", "elink.fcgi", data=data)
    out: dict[str, list[str]] = {}
    for linkset in js.get("linksets") or []:
        ids = linkset.get("ids") or []
        if not ids:
            continue
        key = str(ids[0])
        sra_ids: list[str] = []
        for ldb in linkset.get("linksetdbs") or []:
            if ldb.get("dbto") == "sra":
                sra_ids.extend(str(x) for x in (ldb.get("links") or []))
        out[key] = sra_ids
    return out


def parse_esummary_doc(uid: str, doc: dict) -> dict | None:
    """Parse SRA esummary into structured experiment record."""
    expxml = doc.get("expxml") or doc.get("ExpXml") or ""
    runs_xml = doc.get("runs") or doc.get("Runs") or ""
    m_exp = re.search(r'<Experiment\s+acc="([^"]+)"', expxml)
    if not m_exp:
        return None
    exp = m_exp.group(1)
    title_m = re.search(r"<Title>(.*?)</Title>", expxml)
    title = title_m.group(1) if title_m else exp
    bp_m = re.search(r"<Bioproject>([^<]+)</Bioproject>", expxml)
    bioproject = bp_m.group(1).strip() if bp_m else None
    bs_m = re.search(r"<Biosample>([^<]+)</Biosample>", expxml)
    biosample = bs_m.group(1).strip() if bs_m else None
    study_m = re.search(r'<Study\s+acc="([^"]+)"', expxml)
    study = study_m.group(1) if study_m else None
    plat_m = re.search(r"<Platform[^>]*>([^<]*)</Platform>", expxml)
    instrument_m = re.search(r'instrument_model="([^"]+)"', expxml)
    platform = instrument_m.group(1) if instrument_m else (plat_m.group(1) if plat_m else None)
    stats_m = re.search(
        r'<Statistics[^>]*total_runs="(\d+)"[^>]*total_spots="(\d+)"[^>]*total_bases="(\d+)"',
        expxml,
    )
    spots = int(stats_m.group(2)) if stats_m else None
    bases = int(stats_m.group(3)) if stats_m else None
    run_accs = re.findall(r'<Run\s+acc="([^"]+)"', runs_xml)
    # GSM from title prefix "GSMxxxx:"
    gsm_m = re.match(r"(GSM\d+)\s*:", title)
    gsm = gsm_m.group(1) if gsm_m else None
    lib_m = re.search(r"<LIBRARY_STRATEGY>([^<]+)</LIBRARY_STRATEGY>", expxml)
    return {
        "uid": uid,
        "accession": exp,
        "bioproject": bioproject,
        "biosample": biosample,
        "study": study,
        "title": title,
        "platform": platform,
        "spots": spots,
        "bases": bases,
        "run_accessions": run_accs,
        "gsm": gsm,
        "library_strategy": lib_m.group(1) if lib_m else None,
        # raw_xml omitted
    }


def fetch_esummary_parsed(sra_uids: list[str]) -> list[dict]:
    out: list[dict] = []
    for i in range(0, len(sra_uids), ESUMMARY_BATCH):
        chunk = sra_uids[i : i + ESUMMARY_BATCH]
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
            parsed = parse_esummary_doc(str(uid), doc)
            if parsed:
                out.append(parsed)
    return out


def upsert_experiments(
    conn: sqlite3.Connection,
    dataset_id: str,
    exps: list[dict],
    fetched_at: str,
    stats: dict,
) -> None:
    sample_runs: dict[str, set[str]] = {}
    studies_added: set[str] = set()
    bps_added: set[str] = set()

    for e in exps:
        exp = e["accession"]
        conn.execute(
            """
            INSERT INTO src_sra_experiment(
                accession, bioproject, biosample, title, platform, spots, bases,
                run_accessions, ftp_paths, raw_xml, fetched_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(accession) DO UPDATE SET
                bioproject=excluded.bioproject,
                biosample=excluded.biosample,
                title=excluded.title,
                platform=excluded.platform,
                spots=excluded.spots,
                bases=excluded.bases,
                run_accessions=excluded.run_accessions,
                raw_xml=NULL,
                fetched_at=excluded.fetched_at
            """,
            (
                exp,
                e.get("bioproject"),
                e.get("biosample"),
                e.get("title"),
                e.get("platform"),
                e.get("spots"),
                e.get("bases"),
                json.dumps(e.get("run_accessions") or []),
                None,
                None,  # raw_xml omitted
                fetched_at,
            ),
        )
        stats["experiments_upserted"] = stats.get("experiments_upserted", 0) + 1

        try:
            conn.execute(
                """
                INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
                VALUES (?, 'sra', ?, 'experiment', ?)
                ON CONFLICT(source, accession) DO UPDATE SET
                    dataset_id=excluded.dataset_id,
                    role=excluded.role,
                    external_url=excluded.external_url
                """,
                (dataset_id, exp, f"https://www.ncbi.nlm.nih.gov/sra/{exp}"),
            )
            stats["members_added"] = stats.get("members_added", 0) + 1
        except sqlite3.IntegrityError:
            pass

        study = e.get("study")
        if study and study not in studies_added:
            studies_added.add(study)
            try:
                conn.execute(
                    """
                    INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
                    VALUES (?, 'sra', ?, 'study', ?)
                    ON CONFLICT(source, accession) DO UPDATE SET
                        dataset_id=excluded.dataset_id, role=excluded.role
                    """,
                    (dataset_id, study, f"https://www.ncbi.nlm.nih.gov/sra?term={study}"),
                )
                stats["members_added"] = stats.get("members_added", 0) + 1
            except sqlite3.IntegrityError:
                pass

        bp = e.get("bioproject")
        if bp and bp not in bps_added:
            bps_added.add(bp)
            try:
                conn.execute(
                    """
                    INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
                    VALUES (?, 'sra', ?, 'bioproject', ?)
                    ON CONFLICT(source, accession) DO UPDATE SET
                        dataset_id=excluded.dataset_id, role=excluded.role
                    """,
                    (
                        dataset_id,
                        bp,
                        f"https://www.ncbi.nlm.nih.gov/bioproject/{bp}",
                    ),
                )
                stats["members_added"] = stats.get("members_added", 0) + 1
            except sqlite3.IntegrityError:
                pass

        gsm = e.get("gsm")
        runs = e.get("run_accessions") or []
        if gsm and runs:
            sid = f"geo:{gsm}"
            sample_runs.setdefault(sid, set()).update(runs)

    for sid, runs in sample_runs.items():
        row = conn.execute("SELECT run_ids FROM samples WHERE sample_id=?", (sid,)).fetchone()
        if not row:
            continue
        existing: list[str] = []
        if row[0]:
            try:
                existing = json.loads(row[0]) or []
            except json.JSONDecodeError:
                existing = []
        merged = sorted(set(existing) | runs)
        conn.execute(
            "UPDATE samples SET run_ids=?, updated_at=? WHERE sample_id=?",
            (json.dumps(merged), fetched_at, sid),
        )
        stats["samples_updated"] = stats.get("samples_updated", 0) + 1


def main() -> int:
    conn = connect()
    cur = load_cursor(conn)
    gses = [
        r[0]
        for r in conn.execute(
            "SELECT accession FROM src_geo_series ORDER BY accession"
        ).fetchall()
    ]
    total = len(gses)
    cur["gse_total"] = total
    offset = int(cur.get("offset") or 0)
    log.info("Start geo_sra_link offset=%s total=%s", offset, total)

    already = {
        r[0]
        for r in conn.execute(
            "SELECT DISTINCT dataset_id FROM dataset_members WHERE source='sra' AND role='experiment'"
        )
    }
    log.info("Already have SRA experiment members for %s datasets", len(already))

    batch_since_commit = 0
    i = offset
    try:
        while i < total:
            batch_gses = gses[i : i + ELINK_BATCH]
            uid_map: dict[int, tuple[str, str]] = {}
            to_link: list[int] = []
            for gse in batch_gses:
                ds = f"geo:{gse}"
                if ds in already:
                    cur["gse_processed"] = cur.get("gse_processed", 0) + 1
                    continue
                uid = gds_uid_for_gse(gse)
                if uid is None:
                    cur["gse_processed"] = cur.get("gse_processed", 0) + 1
                    continue
                uid_map[uid] = (gse, ds)
                to_link.append(uid)

            if to_link:
                try:
                    link_map = elink_gds_to_sra(to_link)
                except Exception as e:
                    cur["errors"] = cur.get("errors", 0) + 1
                    log.error("elink batch failed at offset %s: %s", i, e)
                    save_cursor(conn, cur, error=str(e))
                    for _ in batch_gses:
                        cur["gse_processed"] = cur.get("gse_processed", 0) + 1
                    i += len(batch_gses)
                    cur["offset"] = i
                    continue

                sra_to_ds: dict[str, str] = {}
                for uid, (gse, ds) in uid_map.items():
                    sra_ids = link_map.get(str(uid), [])
                    cur["gse_processed"] = cur.get("gse_processed", 0) + 1
                    if not sra_ids:
                        continue
                    cur["gse_with_sra"] = cur.get("gse_with_sra", 0) + 1
                    if len(sra_ids) > MAX_SRA_PER_GSE:
                        cur["truncated_gse"] = cur.get("truncated_gse", 0) + 1
                        sra_ids = sra_ids[:MAX_SRA_PER_GSE]
                    for sid in sra_ids:
                        sra_to_ds[sid] = ds

                if sra_to_ds:
                    try:
                        parsed = fetch_esummary_parsed(list(sra_to_ds.keys()))
                    except Exception as e:
                        cur["errors"] = cur.get("errors", 0) + 1
                        log.error("esummary batch failed at offset %s: %s", i, e)
                        parsed = []

                    by_ds: dict[str, list[dict]] = {}
                    for p in parsed:
                        ds = sra_to_ds.get(p["uid"])
                        if not ds:
                            continue
                        by_ds.setdefault(ds, []).append(p)

                    fetched_at = now_iso()
                    for ds, exps in by_ds.items():
                        try:
                            upsert_experiments(conn, ds, exps, fetched_at, cur)
                            already.add(ds)
                        except Exception as e:
                            cur["errors"] = cur.get("errors", 0) + 1
                            log.error("upsert failed for %s: %s", ds, e)

            i += len(batch_gses)
            cur["offset"] = i
            batch_since_commit += 1
            if batch_since_commit >= COMMIT_EVERY:
                save_cursor(conn, cur)
                batch_since_commit = 0
                log.info(
                    "Progress offset=%s/%s processed=%s with_sra=%s exp=%s samples_upd=%s errors=%s",
                    i,
                    total,
                    cur.get("gse_processed"),
                    cur.get("gse_with_sra"),
                    cur.get("experiments_upserted"),
                    cur.get("samples_updated"),
                    cur.get("errors"),
                )
            elif cur.get("gse_processed", 0) % LOG_EVERY < ELINK_BATCH:
                log.info(
                    "… offset=%s/%s with_sra=%s exp=%s",
                    i,
                    total,
                    cur.get("gse_with_sra"),
                    cur.get("experiments_upserted"),
                )

        cur["offset"] = total
        cur["completed_at"] = now_iso()
        save_cursor(conn, cur)
        log.info("COMPLETED geo_sra_link %s", json.dumps(cur))
        return 0
    except KeyboardInterrupt:
        save_cursor(conn, cur, error="interrupted")
        log.warning("Interrupted at offset %s", cur.get("offset"))
        return 130
    except Exception as e:
        save_cursor(conn, cur, error=str(e))
        log.exception("Fatal: %s", e)
        return 1
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
