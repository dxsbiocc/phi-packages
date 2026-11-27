#!/usr/bin/env python3
"""
GTEx Portal sample-metadata ingest (scan / staging phase).

Writes ONLY to omics_gtex_staging.db — never to omics_v1.db.
Sample / subject annotation only — no expression matrices.
Uses GTEx Portal API v2 (https://gtexportal.org/api/v2/).
Resumable via sync_cursors (source=gtex_sample_scan).
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
STAGING_DB = _DATA_DIR / "omics_gtex_staging.db"
SCHEMA = SCHEMA_PATH
API = "https://gtexportal.org/api/v2"
CURSOR_SOURCE = "gtex_sample_scan"
MIN_INTERVAL = 0.30
MAX_RETRIES = 8
PAGE_SIZE = 500
TZ = timezone(timedelta(hours=8))
SAVE_EVERY = 2  # pages
LOG_EVERY = 2

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "gtex_sample_scan.log"),
    ],
)
log = logging.getLogger("gtex_sample_scan")
_last_req = 0.0


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def throttle() -> None:
    global _last_req
    elapsed = time.monotonic() - _last_req
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_req = time.monotonic()


def gtex_get(path: str, params: dict | None = None) -> Any:
    url = f"{API}/{path.lstrip('/')}"
    last_err: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        throttle()
        try:
            r = http_get(url, params=params, timeout=TIMEOUT)
            if r.status_code == 429 or r.status_code >= 500:
                wait = min(2**attempt, 60)
                log.warning("HTTP %s %s (attempt %s), sleep %ss", r.status_code, path, attempt, wait)
                time.sleep(wait)
                continue
            if r.status_code in (400, 404):
                log.warning("HTTP %s %s body=%s", r.status_code, path, r.text[:300])
                return None
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last_err = e
            wait = min(2**attempt, 60)
            log.warning("Request failed attempt %s %s: %s; sleep %ss", attempt, path, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"GTEx request failed after {MAX_RETRIES}: {last_err}")


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
        CREATE TABLE IF NOT EXISTS src_gtex_subject (
            subject_id             TEXT PRIMARY KEY,
            sex                    TEXT,
            age_bracket            TEXT,
            hardy_scale            TEXT,
            dataset_id             TEXT,
            raw_json               TEXT,
            fetched_at             TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_src_gtex_subject_fetched ON src_gtex_subject(fetched_at);
        CREATE TABLE IF NOT EXISTS src_gtex_sample (
            sample_id              TEXT PRIMARY KEY,  -- e.g. GTEX-1117F-0003-SM-58Q7G
            sample_id_upper        TEXT,
            subject_id             TEXT,
            tissue_sample_id       TEXT,
            aliquot_id             TEXT,
            tissue_site_detail     TEXT,
            tissue_site_detail_id  TEXT,
            uberon_id              TEXT,
            data_type              TEXT,
            sex                    TEXT,
            age_bracket            TEXT,
            hardy_scale            TEXT,
            ischemic_time          REAL,
            ischemic_time_group    TEXT,
            rin                    REAL,
            autolysis_score        TEXT,
            pathology_notes        TEXT,
            pathology_notes_categories TEXT,  -- JSON
            freeze_type            TEXT,
            dataset_id             TEXT,
            raw_json               TEXT,
            fetched_at             TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_src_gtex_sample_tissue ON src_gtex_sample(tissue_site_detail_id);
        CREATE INDEX IF NOT EXISTS idx_src_gtex_sample_subject ON src_gtex_sample(subject_id);
        CREATE INDEX IF NOT EXISTS idx_src_gtex_sample_fetched ON src_gtex_sample(fetched_at);
        CREATE TABLE IF NOT EXISTS src_gtex_tissue (
            tissue_site_detail_id  TEXT PRIMARY KEY,
            tissue_site_detail     TEXT,
            uberon_id              TEXT,
            dataset_id             TEXT,
            color_hex              TEXT,
            raw_json               TEXT,
            fetched_at             TEXT
        );
        CREATE TABLE IF NOT EXISTS gtex_merge_state (
            sample_id TEXT PRIMARY KEY,
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
            "page": 0,
            "total": None,
            "release": None,
            "canonical_dataset": None,
            "subjects_done": False,
            "tissues_done": False,
            "scanned": 0,
            "ingested_samples": 0,
            "ingested_subjects": 0,
            "errors": 0,
            "done": False,
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
        (CURSOR_SOURCE, json.dumps(cur), "page", now_iso(), error),
    )
    conn.commit()


def detect_release() -> tuple[str, str]:
    """Return (api_dataset_id, canonical_id) e.g. ('gtex_v10', 'GTEx:v10')."""
    data = gtex_get("dataset/sample", {"page": 0, "itemsPerPage": 1})
    if not data or not data.get("data"):
        raise RuntimeError("Cannot detect GTEx release from /dataset/sample")
    api_id = data["data"][0].get("datasetId") or "gtex_v10"
    # gtex_v10 -> v10
    m = re.search(r"v(\d+)$", api_id, re.I)
    ver = f"v{m.group(1)}" if m else api_id
    return api_id, f"GTEx:{ver}"


def ensure_dataset(conn: sqlite3.Connection, canonical: str, api_id: str, n_samples: int | None) -> None:
    ts = now_iso()
    existing = conn.execute(
        "SELECT first_seen_at FROM datasets WHERE canonical_id=?", (canonical,)
    ).fetchone()
    first_seen = existing[0] if existing else ts
    title = f"GTEx Analysis {canonical.split(':', 1)[-1]} (Portal {api_id})"
    conn.execute(
        """
        INSERT INTO datasets(
            canonical_id, title, title_norm, organism, tissue, disease_tags,
            assay, platforms, n_samples, access, license, confidence,
            content_hash, first_seen_at, updated_at, probe_status, last_checked
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(canonical_id) DO UPDATE SET
            title=excluded.title,
            title_norm=excluded.title_norm,
            n_samples=COALESCE(excluded.n_samples, datasets.n_samples),
            access=excluded.access,
            confidence=excluded.confidence,
            content_hash=excluded.content_hash,
            updated_at=excluded.updated_at,
            probe_status=excluded.probe_status,
            last_checked=excluded.last_checked
        """,
        (
            canonical,
            title,
            title_norm(title),
            "Homo sapiens",
            "multi-tissue",
            json.dumps([], ensure_ascii=False),
            "RNA-seq",
            json.dumps(["GTEx Portal"], ensure_ascii=False),
            n_samples,
            "public",
            "GTEx Consortium",
            "verified_live",
            content_hash({"release": api_id, "n": n_samples}),
            first_seen,
            ts,
            "ok",
            ts,
        ),
    )
    conn.execute(
        """
        INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
        VALUES (?,?,?,?,?)
        ON CONFLICT(source, accession) DO UPDATE SET
            dataset_id=excluded.dataset_id,
            role=excluded.role,
            external_url=excluded.external_url
        """,
        (
            canonical,
            "gtex",
            api_id,
            "release",
            "https://gtexportal.org/home/",
        ),
    )
    conn.commit()


def ingest_subjects(conn: sqlite3.Connection, cur: dict) -> None:
    if cur.get("subjects_done"):
        return
    page = 0
    total_ingested = 0
    ts = now_iso()
    while True:
        data = gtex_get(
            "dataset/subject",
            {"page": page, "itemsPerPage": PAGE_SIZE, "datasetId": cur.get("release")},
        )
        if data is None:
            raise RuntimeError("subject endpoint returned None")
        items = data.get("data") or []
        paging = data.get("paging_info") or {}
        if page == 0:
            log.info("Subjects total=%s pages=%s", paging.get("totalNumberOfItems"), paging.get("numberOfPages"))
        if not items:
            break
        for sub in items:
            sid = sub.get("subjectId")
            if not sid:
                continue
            conn.execute(
                """
                INSERT INTO src_gtex_subject(
                    subject_id, sex, age_bracket, hardy_scale, dataset_id, raw_json, fetched_at
                ) VALUES (?,?,?,?,?,?,?)
                ON CONFLICT(subject_id) DO UPDATE SET
                    sex=excluded.sex,
                    age_bracket=excluded.age_bracket,
                    hardy_scale=excluded.hardy_scale,
                    dataset_id=excluded.dataset_id,
                    raw_json=excluded.raw_json,
                    fetched_at=excluded.fetched_at
                """,
                (
                    sid,
                    sub.get("sex"),
                    sub.get("ageBracket"),
                    sub.get("hardyScale"),
                    sub.get("datasetId") or cur.get("release"),
                    json.dumps(sub, ensure_ascii=False, default=str),
                    ts,
                ),
            )
            total_ingested += 1
        conn.commit()
        page += 1
        n_pages = paging.get("numberOfPages")
        if n_pages is not None and page >= int(n_pages):
            break
        if len(items) < PAGE_SIZE:
            break
    cur["subjects_done"] = True
    cur["ingested_subjects"] = total_ingested
    save_cursor(conn, cur)
    log.info("Ingested %s GTEx subjects", total_ingested)


def ingest_tissues(conn: sqlite3.Connection, cur: dict) -> None:
    if cur.get("tissues_done"):
        return
    data = gtex_get("dataset/tissueSiteDetail", {"page": 0, "itemsPerPage": 200})
    if data is None:
        log.warning("tissueSiteDetail unavailable; continuing without tissue catalog")
        cur["tissues_done"] = True
        save_cursor(conn, cur)
        return
    items = data.get("data") or []
    ts = now_iso()
    for t in items:
        tid = t.get("tissueSiteDetailId")
        if not tid:
            continue
        conn.execute(
            """
            INSERT INTO src_gtex_tissue(
                tissue_site_detail_id, tissue_site_detail, uberon_id, dataset_id,
                color_hex, raw_json, fetched_at
            ) VALUES (?,?,?,?,?,?,?)
            ON CONFLICT(tissue_site_detail_id) DO UPDATE SET
                tissue_site_detail=excluded.tissue_site_detail,
                uberon_id=excluded.uberon_id,
                dataset_id=excluded.dataset_id,
                color_hex=excluded.color_hex,
                raw_json=excluded.raw_json,
                fetched_at=excluded.fetched_at
            """,
            (
                tid,
                t.get("tissueSiteDetail") or tid.replace("_", " "),
                t.get("uberonId"),
                t.get("datasetId"),
                t.get("colorHex"),
                json.dumps(t, ensure_ascii=False, default=str),
                ts,
            ),
        )
    conn.commit()
    cur["tissues_done"] = True
    save_cursor(conn, cur)
    log.info("Ingested %s GTEx tissue site details", len(items))


def ingest_sample_row(conn: sqlite3.Connection, samp: dict, canonical: str, ts: str) -> int:
    sid = samp.get("sampleId")
    if not sid:
        return 0
    attrs = {
        "sampleId": sid,
        "sampleIdUpper": samp.get("sampleIdUpper"),
        "subjectId": samp.get("subjectId"),
        "tissueSampleId": samp.get("tissueSampleId"),
        "aliquotId": samp.get("aliquotId"),
        "tissueSiteDetail": samp.get("tissueSiteDetail"),
        "tissueSiteDetailId": samp.get("tissueSiteDetailId"),
        "uberonId": samp.get("uberonId"),
        "dataType": samp.get("dataType"),
        "sex": samp.get("sex"),
        "ageBracket": samp.get("ageBracket"),
        "hardyScale": samp.get("hardyScale"),
        "ischemicTime": samp.get("ischemicTime"),
        "ischemicTimeGroup": samp.get("ischemicTimeGroup"),
        "rin": samp.get("rin"),
        "autolysisScore": samp.get("autolysisScore"),
        "pathologyNotes": samp.get("pathologyNotes"),
        "pathologyNotesCategories": samp.get("pathologyNotesCategories"),
        "freezeType": samp.get("freezeType"),
        "datasetId": samp.get("datasetId"),
    }
    path_cats = samp.get("pathologyNotesCategories")
    if not isinstance(path_cats, (dict, list)):
        path_cats_json = None
    else:
        path_cats_json = json.dumps(path_cats, ensure_ascii=False, default=str)

    conn.execute(
        """
        INSERT INTO src_gtex_sample(
            sample_id, sample_id_upper, subject_id, tissue_sample_id, aliquot_id,
            tissue_site_detail, tissue_site_detail_id, uberon_id, data_type,
            sex, age_bracket, hardy_scale, ischemic_time, ischemic_time_group,
            rin, autolysis_score, pathology_notes, pathology_notes_categories,
            freeze_type, dataset_id, raw_json, fetched_at
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(sample_id) DO UPDATE SET
            sample_id_upper=excluded.sample_id_upper,
            subject_id=excluded.subject_id,
            tissue_sample_id=excluded.tissue_sample_id,
            aliquot_id=excluded.aliquot_id,
            tissue_site_detail=excluded.tissue_site_detail,
            tissue_site_detail_id=excluded.tissue_site_detail_id,
            uberon_id=excluded.uberon_id,
            data_type=excluded.data_type,
            sex=excluded.sex,
            age_bracket=excluded.age_bracket,
            hardy_scale=excluded.hardy_scale,
            ischemic_time=excluded.ischemic_time,
            ischemic_time_group=excluded.ischemic_time_group,
            rin=excluded.rin,
            autolysis_score=excluded.autolysis_score,
            pathology_notes=excluded.pathology_notes,
            pathology_notes_categories=excluded.pathology_notes_categories,
            freeze_type=excluded.freeze_type,
            dataset_id=excluded.dataset_id,
            raw_json=excluded.raw_json,
            fetched_at=excluded.fetched_at
        """,
        (
            sid,
            samp.get("sampleIdUpper"),
            samp.get("subjectId"),
            samp.get("tissueSampleId"),
            samp.get("aliquotId"),
            samp.get("tissueSiteDetail"),
            samp.get("tissueSiteDetailId"),
            samp.get("uberonId"),
            samp.get("dataType"),
            samp.get("sex"),
            samp.get("ageBracket"),
            samp.get("hardyScale"),
            samp.get("ischemicTime"),
            samp.get("ischemicTimeGroup"),
            samp.get("rin"),
            str(samp.get("autolysisScore")) if samp.get("autolysisScore") is not None else None,
            samp.get("pathologyNotes"),
            path_cats_json,
            samp.get("freezeType"),
            samp.get("datasetId"),
            json.dumps(samp, ensure_ascii=False, default=str),
            ts,
        ),
    )
    sample_pk = f"GTEx:{sid}"
    tissue = samp.get("tissueSiteDetail") or samp.get("tissueSiteDetailId")
    conn.execute(
        """
        INSERT INTO samples(
            sample_id, dataset_id, title, organism, tissue, disease,
            condition, sex, age, attributes, run_ids, content_hash, updated_at
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(sample_id) DO UPDATE SET
            dataset_id=excluded.dataset_id,
            title=excluded.title,
            organism=excluded.organism,
            tissue=excluded.tissue,
            disease=excluded.disease,
            condition=excluded.condition,
            sex=excluded.sex,
            age=excluded.age,
            attributes=excluded.attributes,
            content_hash=excluded.content_hash,
            updated_at=excluded.updated_at
        """,
        (
            sample_pk,
            canonical,
            sid,
            "Homo sapiens",
            tissue,
            None,  # normal tissue atlas
            samp.get("dataType"),
            samp.get("sex"),
            samp.get("ageBracket"),
            json.dumps(attrs, ensure_ascii=False, default=str),
            None,
            content_hash(attrs),
            ts,
        ),
    )
    return 1


def scan_samples(conn: sqlite3.Connection, cur: dict) -> None:
    if cur.get("done"):
        log.info("Cursor marked done; nothing to scan")
        return
    canonical = cur["canonical_dataset"]
    release = cur["release"]
    page = int(cur.get("page") or 0)
    total = cur.get("total")
    pages_since_save = 0

    while True:
        try:
            data = gtex_get(
                "dataset/sample",
                {"page": page, "itemsPerPage": PAGE_SIZE, "datasetId": release},
            )
        except Exception as e:
            cur["errors"] = cur.get("errors", 0) + 1
            save_cursor(conn, cur, error=str(e)[:500])
            log.error("Fatal page fetch at page=%s: %s", page, e)
            raise

        if data is None:
            cur["errors"] = cur.get("errors", 0) + 1
            save_cursor(conn, cur, error="sample endpoint returned None")
            raise RuntimeError("sample endpoint returned None")

        items = data.get("data") or []
        paging = data.get("paging_info") or {}
        if total is None:
            total = paging.get("totalNumberOfItems")
            cur["total"] = total
            ensure_dataset(conn, canonical, release, total)
            log.info("GTEx samples total=%s pages≈%s release=%s", total, paging.get("numberOfPages"), release)

        if not items:
            log.info("No more samples at page=%s; finishing", page)
            break

        ts = now_iso()
        n = 0
        for samp in items:
            try:
                n += ingest_sample_row(conn, samp, canonical, ts)
            except Exception as e:
                cur["errors"] = cur.get("errors", 0) + 1
                log.warning("Sample ingest error: %s", e)
        conn.commit()

        page += 1
        cur["page"] = page
        cur["scanned"] = cur.get("scanned", 0) + len(items)
        cur["ingested_samples"] = cur.get("ingested_samples", 0) + n
        pages_since_save += 1

        if pages_since_save >= SAVE_EVERY:
            save_cursor(conn, cur)
            pages_since_save = 0

        if page % LOG_EVERY == 0 or (total and cur["scanned"] >= total):
            log.info(
                "progress page=%s scanned=%s/%s ingested=%s errors=%s",
                page,
                cur.get("scanned"),
                total,
                cur.get("ingested_samples"),
                cur.get("errors"),
            )

        n_pages = paging.get("numberOfPages")
        if n_pages is not None and page >= int(n_pages):
            break
        if len(items) < PAGE_SIZE:
            break
        if total is not None and cur.get("scanned", 0) >= total:
            break

    # refresh n_samples
    n_final = conn.execute(
        "SELECT COUNT(*) FROM samples WHERE dataset_id=?", (canonical,)
    ).fetchone()[0]
    conn.execute(
        "UPDATE datasets SET n_samples=?, updated_at=? WHERE canonical_id=?",
        (n_final, now_iso(), canonical),
    )
    conn.commit()

    cur["page"] = page
    cur["done"] = True
    save_cursor(conn, cur)
    log.info(
        "DONE scanned=%s ingested_samples=%s subjects=%s errors=%s",
        cur.get("scanned"),
        cur.get("ingested_samples"),
        cur.get("ingested_subjects"),
        cur.get("errors"),
    )


def main() -> None:
    init_staging()
    conn = connect(STAGING_DB)
    cur = load_cursor(conn)
    log.info("Start cursor=%s", json.dumps(cur))
    try:
        if not cur.get("release") or not cur.get("canonical_dataset"):
            api_id, canonical = detect_release()
            cur["release"] = api_id
            cur["canonical_dataset"] = canonical
            ensure_dataset(conn, canonical, api_id, None)
            save_cursor(conn, cur)
            log.info("Detected GTEx release %s -> %s", api_id, canonical)

        ingest_tissues(conn, cur)
        ingest_subjects(conn, cur)
        scan_samples(conn, cur)
    except Exception as e:
        save_cursor(conn, cur, error=str(e)[:500])
        log.exception("Scan aborted: %s", e)
        conn.close()
        raise
    conn.close()
    log.info("gtex_sample_scan finished cleanly")


if __name__ == "__main__":
    main()
