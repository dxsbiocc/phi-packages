#!/usr/bin/env python3
"""
Remaining PDC programs (non-CPTAC, non-ICPC) sample-metadata ingest via PDC GraphQL.

Writes ONLY to omics_pdc_other_staging.db — never to omics_v1.db.
Sample / case / study metadata + file COUNTS only — no mzML/PSM/protein downloads.
Resumable via sync_cursors (source=pdc_other_scan).

Scope: all PDC programs EXCEPT CPTAC* and ICPC (already staged).
  Covered examples: APOLLO, CBTN, Georgetown, Broad Institute,
  Quantitative digital maps of tissue biopsies.
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

from http_util import TIMEOUT, http_post

from paths import CODE_ROOT as ROOT, SCHEMA_PATH, data_dir as _data_dir, logs_dir as _logs_dir
_DATA_DIR = _data_dir(); _DATA_DIR.mkdir(parents=True, exist_ok=True)
_LOGS_DIR = _logs_dir(); _LOGS_DIR.mkdir(parents=True, exist_ok=True)
STAGING_DB = _DATA_DIR / "omics_pdc_other_staging.db"
SCHEMA = SCHEMA_PATH
API = "https://pdc.cancer.gov/graphql"
CURSOR_SOURCE = "pdc_other_scan"
MIN_INTERVAL = 0.35
MAX_RETRIES = 8
PAGE_SIZE = 50
TZ = timezone(timedelta(hours=8))
SAVE_EVERY_STUDIES = 1

# Already covered by cptac_sample_scan.py / icpc_sample_scan.py — SKIP.
SKIP_PROGRAMS = {
    "Clinical Proteomic Tumor Analysis Consortium",
    "CPTAC-Kids First-Childhood Cancer Data Initiative",
    "CPTAC-Human Cancer Models Initiative",
    "International Cancer Proteogenome Consortium",
}

# Friendly short labels for canonical titles
PROGRAM_SHORT = {
    "Applied Proteogenomics OrganizationaL Learning and Outcomes (APOLLO)": "APOLLO",
    "Pediatric Brain Tumor Atlas - CBTN": "CBTN",
    "Georgetown Proteomics Research Program": "Georgetown",
    "Broad Institute": "Broad",
    "Quantitative digital maps of tissue biopsies": "QDM",
}

LIVER_KW = (
    "liver", "hepat", "hcc", "bile duct", "cholangi", "intrahepatic",
    "gallbladder", "gall bladder", "icc",
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "pdc_other_scan.log"),
    ],
)
log = logging.getLogger("pdc_other_scan")
_last_req = 0.0


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def throttle() -> None:
    global _last_req
    elapsed = time.monotonic() - _last_req
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_req = time.monotonic()


def gql(query: str, timeout: Any = TIMEOUT) -> dict:
    last_err: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        throttle()
        try:
            r = http_post(
                API,
                json={"query": query},
                timeout=timeout,
                headers={"Content-Type": "application/json"},
            )
            if r.status_code == 429 or r.status_code >= 500:
                wait = min(2**attempt, 60)
                log.warning("HTTP %s (attempt %s), sleep %ss", r.status_code, attempt, wait)
                time.sleep(wait)
                continue
            if r.status_code in (400, 404):
                try:
                    body = r.json()
                    errs = body.get("errors") or []
                    msg = errs[0].get("message") if errs else r.text[:300]
                except Exception:
                    msg = r.text[:300]
                log.warning("HTTP %s graphql: %s", r.status_code, msg)
                return {}
            r.raise_for_status()
            data = r.json()
            if data.get("errors"):
                msg = data["errors"][0].get("message", str(data["errors"][0]))
                # treat as soft fail for this query
                log.warning("GraphQL error: %s", msg[:400])
                return {}
            return data.get("data") or {}
        except Exception as e:
            last_err = e
            wait = min(2**attempt, 60)
            log.warning("Request failed attempt %s: %s; sleep %ss", attempt, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"PDC GraphQL failed after {MAX_RETRIES}: {last_err}")


def content_hash(obj: Any) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=False, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:32]


def title_norm(title: str | None) -> str:
    if not title:
        return ""
    return re.sub(r"\s+", " ", title.lower().strip())


def as_list(val: Any) -> list:
    if val is None:
        return []
    if isinstance(val, list):
        return val
    return [val]


def is_liver_blob(*parts: Any) -> bool:
    blob = " ".join(json.dumps(p, ensure_ascii=False, default=str).lower() for p in parts if p is not None)
    return any(k in blob for k in LIVER_KW)


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
        CREATE TABLE IF NOT EXISTS src_pdc_study (
            pdc_study_id           TEXT PRIMARY KEY,
            study_id               TEXT,
            study_submitter_id     TEXT,
            study_name             TEXT,
            program_name           TEXT,
            project_name           TEXT,
            project_submitter_id   TEXT,
            analytical_fraction    TEXT,
            experiment_type        TEXT,
            disease_types          TEXT,   -- JSON
            primary_sites          TEXT,   -- JSON
            cases_count            INTEGER,
            aliquots_count         INTEGER,
            is_liver               INTEGER,
            raw_json               TEXT,
            fetched_at             TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_src_pdc_study_project ON src_pdc_study(project_submitter_id);
        CREATE INDEX IF NOT EXISTS idx_src_pdc_study_liver ON src_pdc_study(is_liver);

        CREATE TABLE IF NOT EXISTS src_pdc_case (
            case_id                TEXT NOT NULL,
            pdc_study_id           TEXT NOT NULL,
            case_submitter_id      TEXT,
            disease_type           TEXT,
            primary_site           TEXT,
            sample_count           INTEGER,
            raw_json               TEXT,
            fetched_at             TEXT,
            PRIMARY KEY (case_id, pdc_study_id)
        );
        CREATE INDEX IF NOT EXISTS idx_src_pdc_case_study ON src_pdc_case(pdc_study_id);
        CREATE INDEX IF NOT EXISTS idx_src_pdc_case_submitter ON src_pdc_case(case_submitter_id);

        CREATE TABLE IF NOT EXISTS src_pdc_sample (
            sample_id              TEXT NOT NULL,
            pdc_study_id           TEXT NOT NULL,
            sample_submitter_id    TEXT,
            case_id                TEXT,
            case_submitter_id      TEXT,
            sample_type            TEXT,
            tissue_type            TEXT,
            gdc_sample_id          TEXT,
            gdc_project_id         TEXT,
            biospecimen_anatomic_site TEXT,
            preservation_method    TEXT,
            disease_type           TEXT,
            primary_site           TEXT,
            analytical_fraction    TEXT,
            experiment_type        TEXT,
            aliquot_count          INTEGER,
            aliquots_json          TEXT,
            raw_json               TEXT,
            fetched_at             TEXT,
            PRIMARY KEY (sample_id, pdc_study_id)
        );
        CREATE INDEX IF NOT EXISTS idx_src_pdc_sample_study ON src_pdc_sample(pdc_study_id);
        CREATE INDEX IF NOT EXISTS idx_src_pdc_sample_case ON src_pdc_sample(case_id);
        CREATE INDEX IF NOT EXISTS idx_src_pdc_sample_gdc ON src_pdc_sample(gdc_sample_id);
        CREATE INDEX IF NOT EXISTS idx_src_pdc_sample_type ON src_pdc_sample(sample_type);

        CREATE TABLE IF NOT EXISTS src_pdc_file_summary (
            pdc_study_id           TEXT PRIMARY KEY,
            by_category            TEXT,   -- JSON [{data_category,file_type,files_count},...]
            total_files            INTEGER,
            fetched_at             TEXT
        );

        CREATE TABLE IF NOT EXISTS pdc_other_merge_state (
            sample_key TEXT PRIMARY KEY,
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
            "studies_done": False,
            "study_index": 0,
            "study_ids": [],
            "current_study_offset": 0,
            "scanned_studies": 0,
            "ingested_studies": 0,
            "ingested_cases": 0,
            "ingested_samples": 0,
            "file_summaries": 0,
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
        (CURSOR_SOURCE, json.dumps(cur), "study_index", now_iso(), error),
    )
    conn.commit()


def fetch_other_studies() -> list[dict]:
    q = """
    {
      allPrograms {
        name
        program_submitter_id
        projects {
          name
          project_submitter_id
          studies {
            pdc_study_id
            study_id
            study_submitter_id
            submitter_id_name
            analytical_fraction
            study_name
            disease_types
            primary_sites
            experiment_type
            cases_count
            aliquots_count
          }
        }
      }
    }
    """
    data = gql(q, timeout=(10, 120))
    programs = data.get("allPrograms") or []
    studies: list[dict] = []
    for p in programs:
        pname = p.get("name") or ""
        # Skip CPTAC family + ICPC (already staged elsewhere)
        if pname in SKIP_PROGRAMS or pname.upper().startswith("CPTAC") or "ICPC" in pname.upper() or "International Cancer Proteogenome" in pname:
            continue
        if not pname:
            continue
        for proj in p.get("projects") or []:
            for s in proj.get("studies") or []:
                if not s.get("pdc_study_id"):
                    continue
                entry = {
                    "program_name": pname,
                    "program_short": PROGRAM_SHORT.get(pname) or pname.split("(")[0].strip()[:40] or "PDC",
                    "project_name": proj.get("name") or proj.get("project_submitter_id"),
                    "project_submitter_id": proj.get("project_submitter_id"),
                    "pdc_study_id": s.get("pdc_study_id"),
                    "study_id": s.get("study_id"),
                    "study_submitter_id": s.get("study_submitter_id"),
                    "submitter_id_name": s.get("submitter_id_name"),
                    "analytical_fraction": s.get("analytical_fraction"),
                    "study_name": s.get("study_name") or s.get("submitter_id_name"),
                    "disease_types": as_list(s.get("disease_types")),
                    "primary_sites": as_list(s.get("primary_sites")),
                    "experiment_type": s.get("experiment_type"),
                    "cases_count": s.get("cases_count"),
                    "aliquots_count": s.get("aliquots_count"),
                }
                studies.append(entry)
    # stable order
    studies.sort(key=lambda x: x["pdc_study_id"])
    # dedupe by pdc_study_id (keep first)
    seen = set()
    uniq = []
    for s in studies:
        if s["pdc_study_id"] in seen:
            continue
        seen.add(s["pdc_study_id"])
        uniq.append(s)
    return uniq


def upsert_study_dataset(conn: sqlite3.Connection, study: dict, ts: str) -> str:
    pid = study["pdc_study_id"]
    canonical = f"PDC:{pid}"
    disease_types = as_list(study.get("disease_types"))
    primary_sites = as_list(study.get("primary_sites"))
    site_join = "; ".join(str(x) for x in primary_sites) if primary_sites else None
    title = study.get("study_name") or study.get("study_submitter_id") or pid
    assay = study.get("analytical_fraction") or "proteomics"
    platforms = [x for x in ["PDC", study.get("program_short") or study.get("program_name"), study.get("experiment_type")] if x]
    liver = 1 if is_liver_blob(disease_types, primary_sites, title) else 0

    slim = {
        k: study.get(k)
        for k in (
            "pdc_study_id",
            "study_id",
            "study_submitter_id",
            "study_name",
            "program_name",
            "project_name",
            "project_submitter_id",
            "analytical_fraction",
            "experiment_type",
            "disease_types",
            "primary_sites",
            "cases_count",
            "aliquots_count",
        )
    }
    conn.execute(
        """
        INSERT INTO src_pdc_study(
            pdc_study_id, study_id, study_submitter_id, study_name,
            program_name, project_name, project_submitter_id,
            analytical_fraction, experiment_type, disease_types, primary_sites,
            cases_count, aliquots_count, is_liver, raw_json, fetched_at
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(pdc_study_id) DO UPDATE SET
            study_id=excluded.study_id,
            study_submitter_id=excluded.study_submitter_id,
            study_name=excluded.study_name,
            program_name=excluded.program_name,
            project_name=excluded.project_name,
            project_submitter_id=excluded.project_submitter_id,
            analytical_fraction=excluded.analytical_fraction,
            experiment_type=excluded.experiment_type,
            disease_types=excluded.disease_types,
            primary_sites=excluded.primary_sites,
            cases_count=excluded.cases_count,
            aliquots_count=excluded.aliquots_count,
            is_liver=excluded.is_liver,
            raw_json=excluded.raw_json,
            fetched_at=excluded.fetched_at
        """,
        (
            pid,
            study.get("study_id"),
            study.get("study_submitter_id"),
            title,
            study.get("program_name"),
            study.get("project_name"),
            study.get("project_submitter_id"),
            study.get("analytical_fraction"),
            study.get("experiment_type"),
            json.dumps(disease_types, ensure_ascii=False),
            json.dumps(primary_sites, ensure_ascii=False),
            study.get("cases_count"),
            study.get("aliquots_count"),
            liver,
            json.dumps(slim, ensure_ascii=False, default=str),
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
            assay=excluded.assay,
            platforms=excluded.platforms,
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
            f"{study.get('program_short') or 'PDC'} {title}",
            title_norm(f"{study.get('program_short') or 'PDC'} {title} {pid}"),
            "Homo sapiens",
            site_join,
            json.dumps(disease_types, ensure_ascii=False),
            assay,
            json.dumps(platforms, ensure_ascii=False),
            study.get("cases_count") or 0,
            "public_metadata",
            "PDC",
            "verified_live",
            content_hash(slim),
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
            "pdc",
            pid,
            "study",
            f"https://pdc.cancer.gov/pdc/study/{pid}",
        ),
    )
    return canonical


def upsert_studies(conn: sqlite3.Connection, cur: dict) -> list[dict]:
    if cur.get("studies_done") and cur.get("study_ids"):
        # reload study metadata from DB for resume
        rows = conn.execute(
            """
            SELECT pdc_study_id, study_id, study_submitter_id, study_name,
                   program_name, project_name, project_submitter_id,
                   analytical_fraction, experiment_type, disease_types,
                   primary_sites, cases_count, aliquots_count
            FROM src_pdc_study
            ORDER BY pdc_study_id
            """
        ).fetchall()
        studies = []
        for r in rows:
            studies.append(
                {
                    "pdc_study_id": r[0],
                    "study_id": r[1],
                    "study_submitter_id": r[2],
                    "study_name": r[3],
                    "program_name": r[4],
                    "program_short": PROGRAM_SHORT.get(r[4] or "") or ((r[4] or "PDC").split("(")[0].strip()[:40] or "PDC"),
                    "project_name": r[5],
                    "project_submitter_id": r[6],
                    "analytical_fraction": r[7],
                    "experiment_type": r[8],
                    "disease_types": json.loads(r[9] or "[]"),
                    "primary_sites": json.loads(r[10] or "[]"),
                    "cases_count": r[11],
                    "aliquots_count": r[12],
                }
            )
        log.info("Resuming with %s PDC-other studies from staging", len(studies))
        return studies

    log.info("Fetching non-CPTAC/non-ICPC study list from PDC allPrograms")
    studies = fetch_other_studies()
    ts = now_iso()
    for s in studies:
        upsert_study_dataset(conn, s, ts)
    conn.commit()
    cur["studies_done"] = True
    cur["study_ids"] = [s["pdc_study_id"] for s in studies]
    cur["ingested_studies"] = len(studies)
    liver_n = sum(1 for s in studies if is_liver_blob(s.get("disease_types"), s.get("primary_sites"), s.get("study_name")))
    log.info(
        "Upserted %s PDC-other studies as datasets (liver-tagged studies=%s)",
        len(studies),
        liver_n,
    )
    save_cursor(conn, cur)
    return studies


def slim_aliquots(aliquots: list) -> list:
    out = []
    for a in aliquots or []:
        if not isinstance(a, dict):
            continue
        out.append(
            {
                "aliquot_id": a.get("aliquot_id"),
                "aliquot_submitter_id": a.get("aliquot_submitter_id"),
                "analyte_type": a.get("analyte_type"),
            }
        )
    return out


def ingest_case_row(
    conn: sqlite3.Connection,
    study: dict,
    case: dict,
    dataset_id: str,
    ts: str,
) -> tuple[int, int]:
    case_id = case.get("case_id")
    if not case_id:
        return 0, 0
    pid = study["pdc_study_id"]
    case_submitter = case.get("case_submitter_id")
    disease_type = case.get("disease_type")
    primary_site = case.get("primary_site")
    samples = case.get("samples") or []

    slim_case = {
        "case_id": case_id,
        "case_submitter_id": case_submitter,
        "disease_type": disease_type,
        "primary_site": primary_site,
        "pdc_study_id": pid,
        "sample_ids": [s.get("sample_id") for s in samples if isinstance(s, dict)],
    }
    conn.execute(
        """
        INSERT INTO src_pdc_case(
            case_id, pdc_study_id, case_submitter_id, disease_type, primary_site,
            sample_count, raw_json, fetched_at
        ) VALUES (?,?,?,?,?,?,?,?)
        ON CONFLICT(case_id, pdc_study_id) DO UPDATE SET
            case_submitter_id=excluded.case_submitter_id,
            disease_type=excluded.disease_type,
            primary_site=excluded.primary_site,
            sample_count=excluded.sample_count,
            raw_json=excluded.raw_json,
            fetched_at=excluded.fetched_at
        """,
        (
            case_id,
            pid,
            case_submitter,
            disease_type,
            primary_site,
            len(samples),
            json.dumps(slim_case, ensure_ascii=False, default=str),
            ts,
        ),
    )

    n_samp = 0
    for samp in samples:
        if not isinstance(samp, dict):
            continue
        sid = samp.get("sample_id")
        if not sid:
            continue
        aliquots = slim_aliquots(samp.get("aliquots") or [])
        samp_submitter = samp.get("sample_submitter_id")
        gdc_sid = samp.get("gdc_sample_id")
        gdc_pid = samp.get("gdc_project_id")
        # normalize null-ish
        if gdc_sid in (None, "null", "None", ""):
            gdc_sid = None
        if gdc_pid in (None, "null", "None", ""):
            gdc_pid = None

        raw_samp = {
            "sample_id": sid,
            "sample_submitter_id": samp_submitter,
            "sample_type": samp.get("sample_type"),
            "tissue_type": samp.get("tissue_type"),
            "gdc_sample_id": gdc_sid,
            "gdc_project_id": gdc_pid,
            "biospecimen_anatomic_site": samp.get("biospecimen_anatomic_site"),
            "preservation_method": samp.get("preservation_method"),
            "composition": samp.get("composition"),
            "aliquots": aliquots,
        }
        conn.execute(
            """
            INSERT INTO src_pdc_sample(
                sample_id, pdc_study_id, sample_submitter_id, case_id, case_submitter_id,
                sample_type, tissue_type, gdc_sample_id, gdc_project_id,
                biospecimen_anatomic_site, preservation_method, disease_type, primary_site,
                analytical_fraction, experiment_type, aliquot_count, aliquots_json,
                raw_json, fetched_at
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            ON CONFLICT(sample_id, pdc_study_id) DO UPDATE SET
                sample_submitter_id=excluded.sample_submitter_id,
                case_id=excluded.case_id,
                case_submitter_id=excluded.case_submitter_id,
                sample_type=excluded.sample_type,
                tissue_type=excluded.tissue_type,
                gdc_sample_id=excluded.gdc_sample_id,
                gdc_project_id=excluded.gdc_project_id,
                biospecimen_anatomic_site=excluded.biospecimen_anatomic_site,
                preservation_method=excluded.preservation_method,
                disease_type=excluded.disease_type,
                primary_site=excluded.primary_site,
                analytical_fraction=excluded.analytical_fraction,
                experiment_type=excluded.experiment_type,
                aliquot_count=excluded.aliquot_count,
                aliquots_json=excluded.aliquots_json,
                raw_json=excluded.raw_json,
                fetched_at=excluded.fetched_at
            """,
            (
                sid,
                pid,
                samp_submitter,
                case_id,
                case_submitter,
                samp.get("sample_type"),
                samp.get("tissue_type"),
                gdc_sid,
                gdc_pid,
                samp.get("biospecimen_anatomic_site"),
                samp.get("preservation_method"),
                disease_type,
                primary_site,
                study.get("analytical_fraction"),
                study.get("experiment_type"),
                len(aliquots),
                json.dumps(aliquots, ensure_ascii=False),
                json.dumps(raw_samp, ensure_ascii=False, default=str),
                ts,
            ),
        )

        sample_pk = f"PDC:{pid}:{sid}"
        attrs = {
            "pdc_study_id": pid,
            "pdc_sample_id": sid,
            "sample_submitter_id": samp_submitter,
            "case_id": case_id,
            "case_submitter_id": case_submitter,
            "sample_type": samp.get("sample_type"),
            "tissue_type": samp.get("tissue_type"),
            "gdc_sample_id": gdc_sid,
            "gdc_project_id": gdc_pid,
            "tcga_submitter_hint": case_submitter if case_submitter and (
                str(case_submitter).startswith("TCGA-") or (gdc_pid and "TCGA" in str(gdc_pid))
            ) else None,
            "biospecimen_anatomic_site": samp.get("biospecimen_anatomic_site"),
            "preservation_method": samp.get("preservation_method"),
            "disease_type": disease_type,
            "primary_site": primary_site,
            "analytical_fraction": study.get("analytical_fraction"),
            "experiment_type": study.get("experiment_type"),
            "program_name": study.get("program_name"),
            "project_submitter_id": study.get("project_submitter_id"),
            "aliquot_count": len(aliquots),
            "aliquot_submitter_ids": [a.get("aliquot_submitter_id") for a in aliquots],
            "is_liver": is_liver_blob(disease_type, primary_site, samp.get("biospecimen_anatomic_site")),
        }
        # drop None attrs for compactness
        attrs = {k: v for k, v in attrs.items() if v is not None}

        tissue = samp.get("tissue_type") or samp.get("biospecimen_anatomic_site") or primary_site
        condition = samp.get("sample_type")
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
                disease_type,
                condition,
                None,
                None,
                json.dumps(attrs, ensure_ascii=False, default=str),
                None,
                content_hash(attrs),
                ts,
            ),
        )
        n_samp += 1
    return 1, n_samp


def fetch_file_summary(conn: sqlite3.Connection, pdc_study_id: str, ts: str) -> None:
    q = """
    {
      filesCountPerStudy(pdc_study_id: "%s") {
        data_category
        file_type
        files_count
      }
    }
    """ % pdc_study_id
    data = gql(q)
    rows = data.get("filesCountPerStudy") or []
    total = sum(int(r.get("files_count") or 0) for r in rows)
    conn.execute(
        """
        INSERT INTO src_pdc_file_summary(pdc_study_id, by_category, total_files, fetched_at)
        VALUES (?,?,?,?)
        ON CONFLICT(pdc_study_id) DO UPDATE SET
            by_category=excluded.by_category,
            total_files=excluded.total_files,
            fetched_at=excluded.fetched_at
        """,
        (pdc_study_id, json.dumps(rows, ensure_ascii=False), total, ts),
    )


def scan_study_biospecimens(
    conn: sqlite3.Connection,
    study: dict,
    cur: dict,
) -> None:
    pid = study["pdc_study_id"]
    dataset_id = f"PDC:{pid}"
    offset = int(cur.get("current_study_offset") or 0)
    total = None

    while True:
        q = """
        {
          paginatedCasesSamplesAliquots(pdc_study_id: "%s" offset: %d limit: %d) {
            total
            casesSamplesAliquots {
              case_id
              case_submitter_id
              disease_type
              primary_site
              samples {
                sample_id
                sample_submitter_id
                sample_type
                tissue_type
                gdc_sample_id
                gdc_project_id
                biospecimen_anatomic_site
                preservation_method
                composition
                aliquots {
                  aliquot_id
                  aliquot_submitter_id
                  analyte_type
                }
              }
            }
          }
        }
        """ % (pid, offset, PAGE_SIZE)
        try:
            data = gql(q, timeout=(10, 120))
        except Exception as e:
            cur["errors"] = cur.get("errors", 0) + 1
            save_cursor(conn, cur, error=f"{pid} offset={offset}: {e}"[:500])
            raise

        block = data.get("paginatedCasesSamplesAliquots") or {}
        if total is None:
            total = block.get("total")
        hits = block.get("casesSamplesAliquots") or []
        if not hits:
            break

        ts = now_iso()
        page_cases = 0
        page_samples = 0
        for case in hits:
            try:
                c, s = ingest_case_row(conn, study, case, dataset_id, ts)
                page_cases += c
                page_samples += s
            except Exception as e:
                cur["errors"] = cur.get("errors", 0) + 1
                log.warning("Case ingest error %s: %s", pid, e)

        conn.commit()
        offset += len(hits)
        cur["current_study_offset"] = offset
        cur["ingested_cases"] = cur.get("ingested_cases", 0) + page_cases
        cur["ingested_samples"] = cur.get("ingested_samples", 0) + page_samples

        if total is not None and offset >= total:
            break
        if len(hits) < PAGE_SIZE:
            break

    # refresh n_samples for this dataset
    n = conn.execute(
        "SELECT COUNT(*) FROM samples WHERE dataset_id=?", (dataset_id,)
    ).fetchone()[0]
    conn.execute(
        "UPDATE datasets SET n_samples=?, updated_at=? WHERE canonical_id=?",
        (n, now_iso(), dataset_id),
    )
    try:
        fetch_file_summary(conn, pid, now_iso())
        cur["file_summaries"] = cur.get("file_summaries", 0) + 1
    except Exception as e:
        cur["errors"] = cur.get("errors", 0) + 1
        log.warning("File summary failed for %s: %s", pid, e)
    conn.commit()
    cur["current_study_offset"] = 0


def scan_all(conn: sqlite3.Connection, cur: dict, studies: list[dict]) -> None:
    if cur.get("done"):
        log.info("Cursor marked done; nothing to scan")
        return
    # map by id
    by_id = {s["pdc_study_id"]: s for s in studies}
    study_ids = cur.get("study_ids") or [s["pdc_study_id"] for s in studies]
    start_idx = int(cur.get("study_index") or 0)

    for i in range(start_idx, len(study_ids)):
        pid = study_ids[i]
        study = by_id.get(pid)
        if not study:
            log.warning("Missing study metadata for %s; skip", pid)
            cur["study_index"] = i + 1
            cur["current_study_offset"] = 0
            continue
        log.info(
            "Scanning study %s/%s %s (%s / %s) cases_count=%s offset=%s",
            i + 1,
            len(study_ids),
            pid,
            study.get("study_name"),
            study.get("analytical_fraction"),
            study.get("cases_count"),
            cur.get("current_study_offset") or 0,
        )
        try:
            scan_study_biospecimens(conn, study, cur)
        except Exception as e:
            save_cursor(conn, cur, error=str(e)[:500])
            log.exception("Aborted on study %s: %s", pid, e)
            raise
        cur["study_index"] = i + 1
        cur["scanned_studies"] = cur.get("scanned_studies", 0) + 1
        cur["current_study_offset"] = 0
        if (i + 1) % SAVE_EVERY_STUDIES == 0:
            save_cursor(conn, cur)
            log.info(
                "progress studies=%s/%s cases=%s samples=%s file_summaries=%s errors=%s",
                i + 1,
                len(study_ids),
                cur.get("ingested_cases"),
                cur.get("ingested_samples"),
                cur.get("file_summaries"),
                cur.get("errors"),
            )

    cur["done"] = True
    save_cursor(conn, cur)
    log.info(
        "DONE studies=%s cases=%s samples=%s file_summaries=%s errors=%s",
        cur.get("scanned_studies"),
        cur.get("ingested_cases"),
        cur.get("ingested_samples"),
        cur.get("file_summaries"),
        cur.get("errors"),
    )


def main() -> None:
    init_staging()
    conn = connect(STAGING_DB)
    cur = load_cursor(conn)
    log.info("Start cursor=%s", json.dumps({k: cur.get(k) for k in cur if k != "study_ids"}))
    try:
        studies = upsert_studies(conn, cur)
        scan_all(conn, cur, studies)
    except Exception as e:
        save_cursor(conn, cur, error=str(e)[:500])
        log.exception("Scan aborted: %s", e)
        conn.close()
        raise
    conn.close()
    log.info("pdc_other_scan finished cleanly")


if __name__ == "__main__":
    main()
