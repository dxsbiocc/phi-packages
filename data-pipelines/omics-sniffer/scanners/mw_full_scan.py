#!/usr/bin/env python3
"""
Metabolomics Workbench full ingest (scan / staging phase).

Writes ONLY to omics_mw_staging.db — never to omics_v1.db.
Resumable via sync_cursors (source=mw_full_scan).
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
from urllib.parse import quote

from http_util import TIMEOUT, http_get

from paths import CODE_ROOT as ROOT, SCHEMA_PATH, data_dir as _data_dir, logs_dir as _logs_dir
_DATA_DIR = _data_dir(); _DATA_DIR.mkdir(parents=True, exist_ok=True)
_LOGS_DIR = _logs_dir(); _LOGS_DIR.mkdir(parents=True, exist_ok=True)
STAGING_DB = _DATA_DIR / "omics_mw_staging.db"
SCHEMA = SCHEMA_PATH
API = "https://www.metabolomicsworkbench.org/rest"
CURSOR_SOURCE = "mw_full_scan"
MIN_INTERVAL = 0.35
MAX_RETRIES = 8
TZ = timezone(timedelta(hours=8))
SAVE_EVERY = 25
LOG_EVERY = 50
BATCH_COMMIT = 10
FETCH_FILES = True
FETCH_ANALYSIS = True
FETCH_DISEASE = True

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "mw_full_scan.log"),
    ],
)
log = logging.getLogger("mw_full_scan")
_last_req = 0.0


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def throttle() -> None:
    global _last_req
    elapsed = time.monotonic() - _last_req
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_req = time.monotonic()


def parse_mw_json(text: str) -> Any:
    """MW sometimes emits invalid JSON escapes in file metadata; recover once."""
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        fixed = re.sub(r'\\(?!["\\/bfnrtu])', r'\\\\', text)
        return json.loads(fixed)


def mw_get(url: str, params: dict | None = None) -> Any:
    last_err: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        throttle()
        try:
            r = http_get(url, params=params, timeout=TIMEOUT)
            if r.status_code == 429 or r.status_code >= 500:
                wait = min(2**attempt, 30)
                log.warning("HTTP %s %s (attempt %s), sleep %ss", r.status_code, url, attempt, wait)
                time.sleep(wait)
                continue
            if r.status_code in (400, 404):
                return None
            r.raise_for_status()
            ctype = r.headers.get("content-type", "")
            text = r.text
            if "json" in ctype or text.lstrip()[:1] in "{[":
                try:
                    return parse_mw_json(text)
                except json.JSONDecodeError as e:
                    # Bad payload is deterministic; do not burn retries.
                    log.warning("JSON parse failed %s: %s", url, e)
                    return None
            return text
        except Exception as e:
            last_err = e
            wait = min(2**attempt, 30)
            log.warning("Request failed attempt %s %s: %s; sleep %ss", attempt, url, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"MW request failed after {MAX_RETRIES}: {last_err}")


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
        CREATE TABLE IF NOT EXISTS src_mw_study (
            accession              TEXT PRIMARY KEY,  -- ST######
            project_id             TEXT,
            title                  TEXT,
            species                TEXT,
            institute              TEXT,
            analysis_type          TEXT,
            number_of_samples      INTEGER,
            submission_date        TEXT,
            release_date           TEXT,
            license                TEXT,
            diseases               TEXT,   -- JSON array
            analyses               TEXT,   -- JSON array of analysis snapshots
            sample_source          TEXT,
            files                  TEXT,   -- JSON compact file list
            file_count             INTEGER,
            raw_data_size_bytes    INTEGER,
            raw_json               TEXT,
            fetched_at             TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_src_mw_fetched_at ON src_mw_study(fetched_at);
        CREATE TABLE IF NOT EXISTS mw_merge_state (
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
            "ingested_new": 0,
            "src_upserted": 0,
            "artifacts_created": 0,
            "errors": 0,
            "files_fetched": 0,
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


def fetch_study_list() -> list[dict]:
    """Unique ST accessions with first-seen project_id from available list."""
    data = mw_get(f"{API}/study/study_title/ST/available")
    if not isinstance(data, dict):
        raise RuntimeError(f"Unexpected MW study list type: {type(data)}")
    by_study: dict[str, dict] = {}
    for v in data.values():
        if not isinstance(v, dict):
            continue
        sid = (v.get("study_id") or "").strip().upper()
        if not sid.startswith("ST"):
            continue
        if sid not in by_study:
            by_study[sid] = {
                "study_id": sid,
                "project_id": (v.get("project_id") or "").strip() or None,
                "analysis_ids": [],
            }
        aid = (v.get("analysis_id") or "").strip()
        if aid and aid not in by_study[sid]["analysis_ids"]:
            by_study[sid]["analysis_ids"].append(aid)

    def sort_key(a: str) -> int:
        m = re.search(r"(\d+)$", a)
        return int(m.group(1)) if m else 0

    ordered = sorted(by_study.keys(), key=sort_key)
    return [by_study[k] for k in ordered]


def as_int(val: Any) -> int | None:
    if val is None or val == "" or val == "-":
        return None
    try:
        return int(float(str(val).replace(",", "")))
    except (TypeError, ValueError):
        return None


def normalize_disease(payload: Any) -> list[str]:
    out: list[str] = []
    if not payload:
        return out
    if isinstance(payload, list):
        for item in payload:
            if isinstance(item, dict):
                name = item.get("Disease") or item.get("disease") or item.get("name")
                if name:
                    out.append(str(name))
            elif item:
                out.append(str(item))
    elif isinstance(payload, dict):
        for item in payload.values():
            if isinstance(item, dict):
                name = item.get("Disease") or item.get("disease") or item.get("name")
                if name:
                    out.append(str(name))
            elif isinstance(item, str) and item:
                out.append(item)
    return list(dict.fromkeys(out))


def compact_files(payload: Any) -> tuple[list[dict], int | None]:
    files: list[dict] = []
    raw_size = None
    if not isinstance(payload, dict):
        return files, raw_size
    top = payload.get("files") or []
    compressed = payload.get("compressed_file_content") or {}
    if isinstance(top, list):
        for name in top:
            if not name:
                continue
            entries = compressed.get(name) if isinstance(compressed, dict) else None
            if isinstance(entries, list) and entries:
                for ent in entries[:200]:
                    if not isinstance(ent, dict):
                        continue
                    files.append(
                        {
                            "archive": name,
                            "name": ent.get("name"),
                            "size": as_int(ent.get("size")),
                        }
                    )
            else:
                files.append({"archive": name, "name": name, "size": None})
    return files, raw_size


def staging_ingest(conn: sqlite3.Connection, entry: dict, cur: dict) -> None:
    accession = entry["study_id"]
    exists = conn.execute(
        "SELECT 1 FROM src_mw_study WHERE accession=? LIMIT 1", (accession,)
    ).fetchone()
    if exists:
        cur["skipped_existing"] = cur.get("skipped_existing", 0) + 1
        return

    fetched_at = now_iso()
    summary = mw_get(f"{API}/study/study_id/{accession}/summary")
    if not isinstance(summary, dict) or not summary.get("study_id"):
        cur["errors"] = cur.get("errors", 0) + 1
        log.warning("no summary for %s", accession)
        return

    title = summary.get("study_title") or accession
    species = summary.get("species")
    institute = summary.get("institute")
    analysis_type = summary.get("analysis_type")
    n_samples = as_int(summary.get("number_of_samples"))
    submission_date = (summary.get("submission_date") or "")[:10] or None
    release_date = (summary.get("release_date") or "")[:10] or None
    license_ = summary.get("license")
    project_id = entry.get("project_id") or summary.get("project_id")

    diseases: list[str] = []
    if FETCH_DISEASE:
        try:
            diseases = normalize_disease(mw_get(f"{API}/study/study_id/{accession}/disease"))
        except Exception as e:
            log.warning("disease fetch failed %s: %s", accession, e)

    sample_source = None
    try:
        src = mw_get(f"{API}/study/study_id/{accession}/source")
        if isinstance(src, dict):
            sample_source = src.get("Sample source") or src.get("sample_source")
    except Exception as e:
        log.warning("source fetch failed %s: %s", accession, e)

    analyses: list[dict] = []
    if FETCH_ANALYSIS:
        try:
            an = mw_get(f"{API}/study/study_id/{accession}/analysis")
            if isinstance(an, dict) and an.get("analysis_id"):
                analyses.append(
                    {
                        "analysis_id": an.get("analysis_id"),
                        "analysis_summary": an.get("analysis_summary"),
                        "analysis_type": an.get("analysis_type"),
                        "ms_instrument_name": an.get("ms_instrument_name"),
                        "ms_instrument_type": an.get("ms_instrument_type"),
                        "chromatography_type": an.get("chromatography_type"),
                        "ion_mode": an.get("ion_mode"),
                    }
                )
            elif isinstance(an, dict):
                for v in an.values():
                    if isinstance(v, dict) and v.get("analysis_id"):
                        analyses.append(
                            {
                                "analysis_id": v.get("analysis_id"),
                                "analysis_summary": v.get("analysis_summary"),
                                "analysis_type": v.get("analysis_type"),
                                "ms_instrument_name": v.get("ms_instrument_name"),
                                "ms_instrument_type": v.get("ms_instrument_type"),
                                "chromatography_type": v.get("chromatography_type"),
                                "ion_mode": v.get("ion_mode"),
                            }
                        )
        except Exception as e:
            log.warning("analysis fetch failed %s: %s", accession, e)

    files_json: list[dict] = []
    raw_size = None
    if FETCH_FILES:
        try:
            files_json, _ = compact_files(mw_get(f"{API}/study/study_id/{accession}/files"))
            cur["files_fetched"] = cur.get("files_fetched", 0) + 1
            rs = mw_get(f"{API}/study/study_id/{accession}/rawdatasize")
            if isinstance(rs, list) and rs:
                for item in rs:
                    if isinstance(item, dict):
                        for v in item.values():
                            raw_size = as_int(v)
                            break
            elif isinstance(rs, dict):
                for v in rs.values():
                    raw_size = as_int(v)
                    break
        except Exception as e:
            log.warning("files fetch failed %s: %s", accession, e)
            cur["errors"] = cur.get("errors", 0) + 1

    raw = {
        "summary": summary,
        "project_id": project_id,
        "analysis_ids": entry.get("analysis_ids") or [],
        "diseases": diseases,
        "sample_source": sample_source,
        "analyses": analyses,
        "n_files": len(files_json),
    }

    conn.execute(
        """
        INSERT INTO src_mw_study(
            accession, project_id, title, species, institute, analysis_type,
            number_of_samples, submission_date, release_date, license, diseases,
            analyses, sample_source, files, file_count, raw_data_size_bytes,
            raw_json, fetched_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(accession) DO UPDATE SET
            project_id=excluded.project_id, title=excluded.title,
            species=excluded.species, institute=excluded.institute,
            analysis_type=excluded.analysis_type,
            number_of_samples=excluded.number_of_samples,
            submission_date=excluded.submission_date,
            release_date=excluded.release_date, license=excluded.license,
            diseases=excluded.diseases, analyses=excluded.analyses,
            sample_source=excluded.sample_source, files=excluded.files,
            file_count=excluded.file_count,
            raw_data_size_bytes=excluded.raw_data_size_bytes,
            raw_json=excluded.raw_json, fetched_at=excluded.fetched_at
        """,
        (
            accession,
            project_id,
            title,
            species,
            institute,
            analysis_type,
            n_samples,
            submission_date,
            release_date,
            license_,
            json.dumps(diseases, ensure_ascii=False),
            json.dumps(analyses, ensure_ascii=False),
            sample_source,
            json.dumps(files_json, ensure_ascii=False),
            len(files_json),
            raw_size,
            json.dumps(raw, ensure_ascii=False, default=str),
            fetched_at,
        ),
    )
    cur["src_upserted"] = cur.get("src_upserted", 0) + 1

    canonical_id = f"mw:{accession}"
    ch = content_hash(
        {"acc": accession, "title": title, "species": species, "n": n_samples, "nf": len(files_json)}
    )
    existing = conn.execute(
        "SELECT first_seen_at FROM datasets WHERE canonical_id=?", (canonical_id,)
    ).fetchone()
    first_seen = existing[0] if existing else fetched_at
    platforms = []
    for a in analyses:
        for k in ("ms_instrument_name", "ms_instrument_type", "chromatography_type"):
            if a.get(k):
                platforms.append(str(a[k]))
    platforms = list(dict.fromkeys(platforms))
    disease_tags = json.dumps(diseases, ensure_ascii=False) if diseases else None

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
            tissue=COALESCE(excluded.tissue, datasets.tissue),
            disease_tags=COALESCE(excluded.disease_tags, datasets.disease_tags),
            assay=COALESCE(excluded.assay, datasets.assay),
            platforms=COALESCE(excluded.platforms, datasets.platforms),
            n_samples=COALESCE(excluded.n_samples, datasets.n_samples),
            license=COALESCE(excluded.license, datasets.license),
            content_hash=excluded.content_hash, updated_at=excluded.updated_at,
            probe_status=excluded.probe_status, last_checked=excluded.last_checked
        """,
        (
            canonical_id,
            title,
            title_norm(title),
            species,
            sample_source,
            disease_tags,
            analysis_type or "metabolomics",
            json.dumps(platforms, ensure_ascii=False) if platforms else None,
            n_samples,
            "public",
            license_,
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

    conn.execute(
        """
        INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
        VALUES (?, 'mw', ?, 'study', ?)
        ON CONFLICT(source, accession) DO UPDATE SET
            dataset_id=excluded.dataset_id, role=excluded.role, external_url=excluded.external_url
        """,
        (
            canonical_id,
            accession,
            summary.get("study_url")
            or f"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID={accession}",
        ),
    )
    if project_id:
        conn.execute(
            """
            INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
            VALUES (?, 'mw_project', ?, 'project', ?)
            ON CONFLICT(source, accession) DO NOTHING
            """,
            (
                canonical_id,
                project_id,
                f"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID={project_id}",
            ),
        )

    conn.execute("DELETE FROM dataset_artifacts WHERE dataset_id=?", (canonical_id,))
    for fi in files_json:
        name = fi.get("name") or fi.get("archive") or "file"
        archive = fi.get("archive") or ""
        kind, role = "supplementary", "study_file"
        lower = str(name).lower()
        if lower.endswith((".raw", ".mzml", ".mzxml", ".wiff", ".d")):
            kind, role = "raw", "raw_ms"
        elif "metadata" in lower or lower.endswith((".xlsx", ".xlsm", ".xls", ".csv", ".tsv")):
            kind, role = "supplementary", "metadata"
        uri = (
            f"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php"
            f"?Mode=Download&StudyID={quote(accession)}&File={quote(str(archive or name))}"
        )
        artifact_id = f"mw:{accession}:{hashlib.sha1(str(name).encode()).hexdigest()[:16]}"
        ach = content_hash({"uri": uri, "name": name, "size": fi.get("size")})
        conn.execute(
            """
            INSERT INTO dataset_artifacts(
                artifact_id, dataset_id, sample_id, kind, role, title, filename,
                format, uri, size_bytes, access, license, source, derived_from,
                probe_status, last_checked, content_hash
            ) VALUES (?, ?, NULL, ?, ?, ?, ?, ?, ?, ?, 'public', ?, 'mw', NULL,
                      'listed', ?, ?)
            ON CONFLICT(artifact_id) DO UPDATE SET
                kind=excluded.kind, role=excluded.role, uri=excluded.uri,
                size_bytes=excluded.size_bytes, last_checked=excluded.last_checked,
                content_hash=excluded.content_hash
            """,
            (
                artifact_id,
                canonical_id,
                kind,
                role,
                name,
                name,
                str(name).rsplit(".", 1)[-1].lower() if "." in str(name) else None,
                uri,
                fi.get("size"),
                license_,
                fetched_at,
                ach,
            ),
        )
        cur["artifacts_created"] = cur.get("artifacts_created", 0) + 1


def main() -> int:
    log.info("Starting Metabolomics Workbench SCAN-ONLY (staging=%s)", STAGING_DB)
    init_staging()
    st = connect(STAGING_DB, busy_ms=5000)
    cur = load_cursor(st)

    try:
        log.info("Fetching MW study accession list...")
        studies = fetch_study_list()
        total = len(studies)
        cur["total"] = total
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
            entry = studies[idx]
            try:
                staging_ingest(st, entry, cur)
                pending_commit += 1
            except Exception as e:
                cur["errors"] = cur.get("errors", 0) + 1
                log.exception("Failed staging %s: %s", entry.get("study_id"), e)
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
                    "progress index=%s/%s scanned=%s new=%s skipped=%s src=%s artifacts=%s files=%s errors=%s",
                    idx,
                    total,
                    cur.get("scanned"),
                    cur.get("ingested_new"),
                    cur.get("skipped_existing"),
                    cur.get("src_upserted"),
                    cur.get("artifacts_created"),
                    cur.get("files_fetched"),
                    cur.get("errors"),
                )

        if pending_commit:
            st.commit()
        cur["index"] = total
        cur["completed_at"] = now_iso()
        save_cursor(st, cur)
        log.info("Finished MW scan: %s", json.dumps(cur))
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
