#!/usr/bin/env python3
"""
GDC / TCGA sample-metadata ingest (scan / staging phase).

Writes ONLY to omics_gdc_tcga_staging.db — never to omics_v1.db.
Sample / case clinical metadata only — no BAM/FASTQ/expression download.
Resumable via sync_cursors (source=gdc_tcga_sample_scan).
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

from http_util import TIMEOUT, http_get, http_post

from paths import CODE_ROOT as ROOT, SCHEMA_PATH, data_dir as _data_dir, logs_dir as _logs_dir
_DATA_DIR = _data_dir(); _DATA_DIR.mkdir(parents=True, exist_ok=True)
_LOGS_DIR = _logs_dir(); _LOGS_DIR.mkdir(parents=True, exist_ok=True)
STAGING_DB = _DATA_DIR / "omics_gdc_tcga_staging.db"
SCHEMA = SCHEMA_PATH
API = "https://api.gdc.cancer.gov"
CURSOR_SOURCE = "gdc_tcga_sample_scan"
MIN_INTERVAL = 0.35
MAX_RETRIES = 8
PAGE_SIZE = 50
TZ = timezone(timedelta(hours=8))
SAVE_EVERY = 5  # pages
LOG_EVERY = 5

CASE_FIELDS = [
    "case_id",
    "submitter_id",
    "project.project_id",
    "project.program.name",
    "primary_site",
    "disease_type",
    "demographic.gender",
    "demographic.race",
    "demographic.ethnicity",
    "demographic.year_of_birth",
    "demographic.vital_status",
    "demographic.days_to_death",
    "diagnoses.primary_diagnosis",
    "diagnoses.ajcc_pathologic_stage",
    "diagnoses.tumor_stage",
    "diagnoses.age_at_diagnosis",
    "diagnoses.vital_status",
    "diagnoses.days_to_last_follow_up",
    "diagnoses.morphology",
    "diagnoses.site_of_resection_or_biopsy",
    "samples.sample_id",
    "samples.submitter_id",
    "samples.sample_type",
    "samples.sample_type_id",
    "samples.tissue_type",
    "samples.tumor_code",
    "samples.tumor_code_id",
    "samples.is_ffpe",
    "samples.oct_embedded",
    "samples.shortest_dimension",
    "samples.longest_dimension",
    "samples.intermediate_dimension",
    "samples.preservation_method",
    "samples.days_to_collection",
    "samples.initial_weight",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "gdc_tcga_sample_scan.log"),
    ],
)
log = logging.getLogger("gdc_tcga_sample_scan")
_last_req = 0.0


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def throttle() -> None:
    global _last_req
    elapsed = time.monotonic() - _last_req
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_req = time.monotonic()


def gdc_request(method: str, path: str, payload: dict | None = None) -> dict:
    url = f"{API}{path}"
    last_err: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        throttle()
        try:
            if method == "GET":
                r = http_get(url, params=payload, timeout=TIMEOUT)
            else:
                r = http_post(url, json=payload, timeout=TIMEOUT)
            if r.status_code == 429 or r.status_code >= 500:
                wait = min(2**attempt, 60)
                log.warning("HTTP %s %s (attempt %s), sleep %ss", r.status_code, path, attempt, wait)
                time.sleep(wait)
                continue
            if r.status_code in (400, 404):
                log.warning("HTTP %s %s body=%s", r.status_code, path, r.text[:300])
                return {}
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last_err = e
            wait = min(2**attempt, 60)
            log.warning("Request failed attempt %s %s: %s; sleep %ss", attempt, path, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"GDC request failed after {MAX_RETRIES}: {last_err}")


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
        CREATE TABLE IF NOT EXISTS src_gdc_project (
            project_id             TEXT PRIMARY KEY,
            name                   TEXT,
            program_name           TEXT,
            primary_site           TEXT,   -- JSON
            disease_type           TEXT,   -- JSON
            case_count             INTEGER,
            file_count             INTEGER,
            raw_json               TEXT,
            fetched_at             TEXT
        );
        CREATE TABLE IF NOT EXISTS src_gdc_case (
            case_id                TEXT PRIMARY KEY,
            submitter_id           TEXT,
            project_id             TEXT,
            primary_site           TEXT,
            disease_type           TEXT,
            gender                 TEXT,
            race                   TEXT,
            ethnicity              TEXT,
            vital_status           TEXT,
            year_of_birth          INTEGER,
            days_to_death          REAL,
            diagnoses_json         TEXT,
            sample_count           INTEGER,
            raw_json               TEXT,
            fetched_at             TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_src_gdc_case_project ON src_gdc_case(project_id);
        CREATE INDEX IF NOT EXISTS idx_src_gdc_case_fetched ON src_gdc_case(fetched_at);
        CREATE TABLE IF NOT EXISTS src_gdc_sample (
            sample_id              TEXT PRIMARY KEY,  -- GDC UUID
            submitter_id           TEXT,
            case_id                TEXT,
            case_submitter_id      TEXT,
            project_id             TEXT,
            sample_type            TEXT,
            sample_type_id         TEXT,
            tissue_type            TEXT,
            tumor_code             TEXT,
            tumor_code_id          TEXT,
            is_ffpe                TEXT,
            oct_embedded           TEXT,
            preservation_method    TEXT,
            days_to_collection     REAL,
            initial_weight         REAL,
            shortest_dimension     REAL,
            longest_dimension      REAL,
            intermediate_dimension REAL,
            raw_json               TEXT,
            fetched_at             TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_src_gdc_sample_project ON src_gdc_sample(project_id);
        CREATE INDEX IF NOT EXISTS idx_src_gdc_sample_case ON src_gdc_sample(case_id);
        CREATE INDEX IF NOT EXISTS idx_src_gdc_sample_fetched ON src_gdc_sample(fetched_at);
        CREATE TABLE IF NOT EXISTS src_gdc_file_summary (
            project_id             TEXT PRIMARY KEY,
            by_data_type           TEXT,   -- JSON {type: count}
            by_data_format         TEXT,   -- JSON {format: count}
            total_files            INTEGER,
            fetched_at             TEXT
        );
        CREATE TABLE IF NOT EXISTS gdc_tcga_merge_state (
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
            "from": 0,
            "total": None,
            "projects_done": False,
            "file_summary_done": False,
            "scanned_cases": 0,
            "ingested_cases": 0,
            "ingested_samples": 0,
            "datasets_upserted": 0,
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
        (CURSOR_SOURCE, json.dumps(cur), "from_offset", now_iso(), error),
    )
    conn.commit()


def as_list(val: Any) -> list:
    if val is None:
        return []
    if isinstance(val, list):
        return val
    return [val]


def slim_case_json(case: dict) -> dict:
    """Drop unused nested blobs; keep clinical + samples."""
    keep = {
        "case_id": case.get("case_id") or case.get("id"),
        "submitter_id": case.get("submitter_id"),
        "primary_site": case.get("primary_site"),
        "disease_type": case.get("disease_type"),
        "project": case.get("project"),
        "demographic": case.get("demographic"),
        "diagnoses": case.get("diagnoses"),
        "samples": case.get("samples"),
    }
    return keep


def upsert_projects(conn: sqlite3.Connection, cur: dict) -> None:
    if cur.get("projects_done"):
        return
    log.info("Fetching TCGA project list from GDC")
    payload = {
        "filters": {"op": "=", "content": {"field": "program.name", "value": "TCGA"}},
        "fields": "project_id,name,primary_site,disease_type,summary.case_count,summary.file_count,program.name",
        "format": "json",
        "size": 100,
    }
    data = gdc_request("POST", "/projects", payload)
    hits = (data.get("data") or {}).get("hits") or []
    ts = now_iso()
    for h in hits:
        pid = h.get("project_id")
        if not pid:
            continue
        summary = h.get("summary") or {}
        primary_site = as_list(h.get("primary_site"))
        disease_type = as_list(h.get("disease_type"))
        case_count = summary.get("case_count")
        file_count = summary.get("file_count")
        name = h.get("name") or pid
        program = (h.get("program") or {}).get("name") or "TCGA"
        site_join = "; ".join(str(x) for x in primary_site) if primary_site else None
        disease_join = "; ".join(str(x) for x in disease_type) if disease_type else None
        canonical = f"TCGA:{pid}"
        conn.execute(
            """
            INSERT INTO src_gdc_project(
                project_id, name, program_name, primary_site, disease_type,
                case_count, file_count, raw_json, fetched_at
            ) VALUES (?,?,?,?,?,?,?,?,?)
            ON CONFLICT(project_id) DO UPDATE SET
                name=excluded.name,
                program_name=excluded.program_name,
                primary_site=excluded.primary_site,
                disease_type=excluded.disease_type,
                case_count=excluded.case_count,
                file_count=excluded.file_count,
                raw_json=excluded.raw_json,
                fetched_at=excluded.fetched_at
            """,
            (
                pid,
                name,
                program,
                json.dumps(primary_site, ensure_ascii=False),
                json.dumps(disease_type, ensure_ascii=False),
                case_count,
                file_count,
                json.dumps(h, ensure_ascii=False, default=str),
                ts,
            ),
        )
        existing = conn.execute(
            "SELECT first_seen_at FROM datasets WHERE canonical_id=?", (canonical,)
        ).fetchone()
        first_seen = existing[0] if existing else ts
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
                tissue=excluded.tissue,
                disease_tags=excluded.disease_tags,
                n_samples=excluded.n_samples,
                access=excluded.access,
                confidence=excluded.confidence,
                content_hash=excluded.content_hash,
                updated_at=excluded.updated_at,
                probe_status=excluded.probe_status,
                last_checked=excluded.last_checked
            """,
            (
                canonical,
                f"TCGA {name} ({pid})",
                title_norm(f"TCGA {name} {pid}"),
                "Homo sapiens",
                site_join,
                json.dumps(disease_type, ensure_ascii=False),
                "multiomics",
                json.dumps(["GDC", "TCGA"], ensure_ascii=False),
                case_count,  # updated later with actual sample counts
                "public_metadata",
                "TCGA/GDC",
                "verified_live",
                content_hash(h),
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
                "gdc",
                pid,
                "project",
                f"https://portal.gdc.cancer.gov/projects/{pid}",
            ),
        )
        cur["datasets_upserted"] = cur.get("datasets_upserted", 0) + 1
    conn.commit()
    cur["projects_done"] = True
    cur["total_projects"] = len(hits)
    log.info("Upserted %s TCGA projects as datasets", len(hits))
    save_cursor(conn, cur)


def upsert_file_summaries(conn: sqlite3.Connection, cur: dict) -> None:
    """Light per-project file counts by data_type / data_format — no downloads."""
    if cur.get("file_summary_done"):
        return
    projects = [
        r[0]
        for r in conn.execute("SELECT project_id FROM src_gdc_project ORDER BY project_id").fetchall()
    ]
    ts = now_iso()
    for pid in projects:
        try:
            # aggregation via facets
            payload = {
                "filters": {
                    "op": "and",
                    "content": [
                        {"op": "=", "content": {"field": "cases.project.project_id", "value": pid}},
                    ],
                },
                "size": 0,
                "facets": "data_type,data_format",
            }
            data = gdc_request("POST", "/files", payload)
            aggs = (data.get("data") or {}).get("aggregations") or {}
            by_type: dict[str, int] = {}
            by_fmt: dict[str, int] = {}
            for bucket in (aggs.get("data_type") or {}).get("buckets") or []:
                by_type[str(bucket.get("key"))] = int(bucket.get("doc_count") or 0)
            for bucket in (aggs.get("data_format") or {}).get("buckets") or []:
                by_fmt[str(bucket.get("key"))] = int(bucket.get("doc_count") or 0)
            total = sum(by_type.values()) or sum(by_fmt.values())
            conn.execute(
                """
                INSERT INTO src_gdc_file_summary(
                    project_id, by_data_type, by_data_format, total_files, fetched_at
                ) VALUES (?,?,?,?,?)
                ON CONFLICT(project_id) DO UPDATE SET
                    by_data_type=excluded.by_data_type,
                    by_data_format=excluded.by_data_format,
                    total_files=excluded.total_files,
                    fetched_at=excluded.fetched_at
                """,
                (
                    pid,
                    json.dumps(by_type, ensure_ascii=False),
                    json.dumps(by_fmt, ensure_ascii=False),
                    total,
                    ts,
                ),
            )
            conn.commit()
            log.info("File summary %s: total=%s types=%s", pid, total, len(by_type))
        except Exception as e:
            cur["errors"] = cur.get("errors", 0) + 1
            log.warning("File summary failed for %s: %s", pid, e)
            save_cursor(conn, cur, error=str(e)[:500])
    cur["file_summary_done"] = True
    save_cursor(conn, cur)


def ingest_case(conn: sqlite3.Connection, case: dict, ts: str) -> tuple[int, int]:
    """Returns (cases_upserted, samples_upserted)."""
    case_id = case.get("case_id") or case.get("id")
    if not case_id:
        return 0, 0
    project = case.get("project") or {}
    project_id = project.get("project_id") or "UNKNOWN"
    dataset_id = f"TCGA:{project_id}"
    submitter_id = case.get("submitter_id")
    demo = case.get("demographic") or {}
    diagnoses = case.get("diagnoses") or []
    samples = case.get("samples") or []
    primary_site = case.get("primary_site")
    if isinstance(primary_site, list):
        primary_site = "; ".join(str(x) for x in primary_site)
    disease_type = case.get("disease_type")
    if isinstance(disease_type, list):
        disease_type = "; ".join(str(x) for x in disease_type)

    # primary diagnosis for disease field
    primary_dx = None
    stage = None
    age_at_dx = None
    if diagnoses:
        d0 = diagnoses[0] if isinstance(diagnoses[0], dict) else {}
        primary_dx = d0.get("primary_diagnosis")
        stage = d0.get("ajcc_pathologic_stage") or d0.get("tumor_stage")
        age_at_dx = d0.get("age_at_diagnosis")

    slim = slim_case_json(case)
    conn.execute(
        """
        INSERT INTO src_gdc_case(
            case_id, submitter_id, project_id, primary_site, disease_type,
            gender, race, ethnicity, vital_status, year_of_birth, days_to_death,
            diagnoses_json, sample_count, raw_json, fetched_at
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(case_id) DO UPDATE SET
            submitter_id=excluded.submitter_id,
            project_id=excluded.project_id,
            primary_site=excluded.primary_site,
            disease_type=excluded.disease_type,
            gender=excluded.gender,
            race=excluded.race,
            ethnicity=excluded.ethnicity,
            vital_status=excluded.vital_status,
            year_of_birth=excluded.year_of_birth,
            days_to_death=excluded.days_to_death,
            diagnoses_json=excluded.diagnoses_json,
            sample_count=excluded.sample_count,
            raw_json=excluded.raw_json,
            fetched_at=excluded.fetched_at
        """,
        (
            case_id,
            submitter_id,
            project_id,
            primary_site,
            disease_type,
            demo.get("gender") or demo.get("sex"),
            demo.get("race"),
            demo.get("ethnicity"),
            demo.get("vital_status"),
            demo.get("year_of_birth"),
            demo.get("days_to_death"),
            json.dumps(diagnoses, ensure_ascii=False, default=str),
            len(samples),
            json.dumps(slim, ensure_ascii=False, default=str),
            ts,
        ),
    )

    # ensure dataset exists (in case projects step missed)
    existing_ds = conn.execute(
        "SELECT first_seen_at FROM datasets WHERE canonical_id=?", (dataset_id,)
    ).fetchone()
    if not existing_ds:
        conn.execute(
            """
            INSERT INTO datasets(
                canonical_id, title, title_norm, organism, tissue, disease_tags,
                assay, platforms, n_samples, access, license, confidence,
                content_hash, first_seen_at, updated_at, probe_status, last_checked
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            """,
            (
                dataset_id,
                f"TCGA {project_id}",
                title_norm(f"TCGA {project_id}"),
                "Homo sapiens",
                primary_site,
                json.dumps([disease_type] if disease_type else [], ensure_ascii=False),
                "multiomics",
                json.dumps(["GDC", "TCGA"], ensure_ascii=False),
                0,
                "public_metadata",
                "TCGA/GDC",
                "verified_live",
                content_hash({"project_id": project_id}),
                ts,
                ts,
                "ok",
                ts,
            ),
        )

    n_samp = 0
    gender = demo.get("gender") or demo.get("sex")
    age_str = None
    if age_at_dx is not None:
        try:
            age_str = f"{float(age_at_dx) / 365.25:.1f}y"
        except (TypeError, ValueError):
            age_str = str(age_at_dx)

    for samp in samples:
        if not isinstance(samp, dict):
            continue
        sid = samp.get("sample_id")
        if not sid:
            continue
        samp_submitter = samp.get("submitter_id")
        sample_pk = f"GDC:{sid}"
        attrs = {
            "gdc_sample_id": sid,
            "submitter_id": samp_submitter,
            "case_id": case_id,
            "case_submitter_id": submitter_id,
            "project_id": project_id,
            "sample_type": samp.get("sample_type"),
            "sample_type_id": samp.get("sample_type_id"),
            "tissue_type": samp.get("tissue_type"),
            "tumor_code": samp.get("tumor_code"),
            "tumor_code_id": samp.get("tumor_code_id"),
            "is_ffpe": samp.get("is_ffpe"),
            "oct_embedded": samp.get("oct_embedded"),
            "preservation_method": samp.get("preservation_method"),
            "days_to_collection": samp.get("days_to_collection"),
            "initial_weight": samp.get("initial_weight"),
            "shortest_dimension": samp.get("shortest_dimension"),
            "longest_dimension": samp.get("longest_dimension"),
            "intermediate_dimension": samp.get("intermediate_dimension"),
            "primary_diagnosis": primary_dx,
            "ajcc_pathologic_stage": stage,
            "age_at_diagnosis_days": age_at_dx,
            "race": demo.get("race"),
            "ethnicity": demo.get("ethnicity"),
            "vital_status": demo.get("vital_status"),
            "days_to_death": demo.get("days_to_death"),
        }
        conn.execute(
            """
            INSERT INTO src_gdc_sample(
                sample_id, submitter_id, case_id, case_submitter_id, project_id,
                sample_type, sample_type_id, tissue_type, tumor_code, tumor_code_id,
                is_ffpe, oct_embedded, preservation_method, days_to_collection,
                initial_weight, shortest_dimension, longest_dimension,
                intermediate_dimension, raw_json, fetched_at
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            ON CONFLICT(sample_id) DO UPDATE SET
                submitter_id=excluded.submitter_id,
                case_id=excluded.case_id,
                case_submitter_id=excluded.case_submitter_id,
                project_id=excluded.project_id,
                sample_type=excluded.sample_type,
                sample_type_id=excluded.sample_type_id,
                tissue_type=excluded.tissue_type,
                tumor_code=excluded.tumor_code,
                tumor_code_id=excluded.tumor_code_id,
                is_ffpe=excluded.is_ffpe,
                oct_embedded=excluded.oct_embedded,
                preservation_method=excluded.preservation_method,
                days_to_collection=excluded.days_to_collection,
                initial_weight=excluded.initial_weight,
                shortest_dimension=excluded.shortest_dimension,
                longest_dimension=excluded.longest_dimension,
                intermediate_dimension=excluded.intermediate_dimension,
                raw_json=excluded.raw_json,
                fetched_at=excluded.fetched_at
            """,
            (
                sid,
                samp_submitter,
                case_id,
                submitter_id,
                project_id,
                samp.get("sample_type"),
                str(samp.get("sample_type_id")) if samp.get("sample_type_id") is not None else None,
                samp.get("tissue_type"),
                samp.get("tumor_code"),
                str(samp.get("tumor_code_id")) if samp.get("tumor_code_id") is not None else None,
                str(samp.get("is_ffpe")) if samp.get("is_ffpe") is not None else None,
                samp.get("oct_embedded"),
                samp.get("preservation_method"),
                samp.get("days_to_collection"),
                samp.get("initial_weight"),
                samp.get("shortest_dimension"),
                samp.get("longest_dimension"),
                samp.get("intermediate_dimension"),
                json.dumps(samp, ensure_ascii=False, default=str),
                ts,
            ),
        )
        tissue = samp.get("tissue_type") or primary_site
        condition = samp.get("sample_type")
        disease = primary_dx or disease_type
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
                dataset_id,
                samp_submitter or sid,
                "Homo sapiens",
                tissue,
                disease,
                condition,
                gender,
                age_str,
                json.dumps(attrs, ensure_ascii=False, default=str),
                None,
                content_hash(attrs),
                ts,
            ),
        )
        n_samp += 1
    return 1, n_samp


def refresh_dataset_sample_counts(conn: sqlite3.Connection) -> None:
    rows = conn.execute(
        """
        SELECT dataset_id, COUNT(*) FROM samples
        WHERE dataset_id LIKE 'TCGA:%'
        GROUP BY dataset_id
        """
    ).fetchall()
    ts = now_iso()
    for did, n in rows:
        conn.execute(
            "UPDATE datasets SET n_samples=?, updated_at=? WHERE canonical_id=?",
            (n, ts, did),
        )
    conn.commit()


def scan_cases(conn: sqlite3.Connection, cur: dict) -> None:
    if cur.get("done"):
        log.info("Cursor marked done; nothing to scan")
        return
    offset = int(cur.get("from") or 0)
    total = cur.get("total")
    pages_since_save = 0

    while True:
        payload = {
            "filters": {
                "op": "=",
                "content": {"field": "project.program.name", "value": "TCGA"},
            },
            "fields": ",".join(CASE_FIELDS),
            "format": "json",
            "size": PAGE_SIZE,
            "from": offset,
        }
        try:
            data = gdc_request("POST", "/cases", payload)
        except Exception as e:
            cur["errors"] = cur.get("errors", 0) + 1
            save_cursor(conn, cur, error=str(e)[:500])
            log.error("Fatal page fetch at from=%s: %s", offset, e)
            raise

        pdata = data.get("data") or {}
        pagination = pdata.get("pagination") or {}
        hits = pdata.get("hits") or []
        if total is None:
            total = pagination.get("total")
            cur["total"] = total
            log.info("TCGA cases total=%s", total)

        if not hits:
            log.info("No more hits at from=%s; finishing", offset)
            break

        ts = now_iso()
        page_cases = 0
        page_samples = 0
        for case in hits:
            try:
                c, s = ingest_case(conn, case, ts)
                page_cases += c
                page_samples += s
            except Exception as e:
                cur["errors"] = cur.get("errors", 0) + 1
                log.warning("Case ingest error: %s", e)

        conn.commit()
        offset += len(hits)
        cur["from"] = offset
        cur["scanned_cases"] = cur.get("scanned_cases", 0) + len(hits)
        cur["ingested_cases"] = cur.get("ingested_cases", 0) + page_cases
        cur["ingested_samples"] = cur.get("ingested_samples", 0) + page_samples
        pages_since_save += 1

        if pages_since_save >= SAVE_EVERY or (total and offset >= total):
            save_cursor(conn, cur)
            pages_since_save = 0

        if cur.get("scanned_cases", 0) % (PAGE_SIZE * LOG_EVERY) < PAGE_SIZE or (
            total and offset >= total
        ):
            log.info(
                "progress from=%s/%s scanned_cases=%s samples=%s errors=%s",
                offset,
                total,
                cur.get("scanned_cases"),
                cur.get("ingested_samples"),
                cur.get("errors"),
            )

        if total is not None and offset >= total:
            break
        if len(hits) < PAGE_SIZE:
            break

    refresh_dataset_sample_counts(conn)
    cur["from"] = offset
    cur["done"] = True
    save_cursor(conn, cur)
    log.info(
        "DONE scanned_cases=%s ingested_cases=%s ingested_samples=%s errors=%s",
        cur.get("scanned_cases"),
        cur.get("ingested_cases"),
        cur.get("ingested_samples"),
        cur.get("errors"),
    )


def main() -> None:
    init_staging()
    conn = connect(STAGING_DB)
    cur = load_cursor(conn)
    log.info("Start cursor=%s", json.dumps(cur))
    try:
        upsert_projects(conn, cur)
        # File summaries are optional and can be slow; do after samples if not done
        scan_cases(conn, cur)
        try:
            upsert_file_summaries(conn, cur)
        except Exception as e:
            log.warning("File summary phase incomplete: %s", e)
            save_cursor(conn, cur, error=f"file_summary: {e}"[:500])
    except Exception as e:
        save_cursor(conn, cur, error=str(e)[:500])
        log.exception("Scan aborted: %s", e)
        conn.close()
        raise
    conn.close()
    log.info("gdc_tcga_sample_scan finished cleanly")


if __name__ == "__main__":
    main()
