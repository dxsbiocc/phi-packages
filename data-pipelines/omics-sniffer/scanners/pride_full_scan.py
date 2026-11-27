#!/usr/bin/env python3
"""
PRIDE Archive full ingest (scan / staging phase).

Writes ONLY to omics_pride_staging.db — never to omics_v1.db —
so SRA orphan scan and ArrayExpress staging can keep the main DB / AE staging
without lock contention.

Resumable via sync_cursors (source=pride_full_scan).
Merge into main is a separate offline step (not implemented here).
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

import requests

from http_util import TIMEOUT, http_get

from paths import CODE_ROOT as ROOT, SCHEMA_PATH, data_dir as _data_dir, logs_dir as _logs_dir
_DATA_DIR = _data_dir(); _DATA_DIR.mkdir(parents=True, exist_ok=True)
_LOGS_DIR = _logs_dir(); _LOGS_DIR.mkdir(parents=True, exist_ok=True)
STAGING_DB = _DATA_DIR / "omics_pride_staging.db"
SCHEMA = SCHEMA_PATH
API = "https://www.ebi.ac.uk/pride/ws/archive/v2"
CURSOR_SOURCE = "pride_full_scan"
PAGE_SIZE = 100
MIN_INTERVAL = 0.28
MAX_RETRIES = 8
TZ = timezone(timedelta(hours=8))
FETCH_FILES = True
MAX_FILES_PER_PROJECT = 300
FILE_PAGE_SIZE = 100
# Drop huge protocol blobs from stored raw_json
RAW_DROP_KEYS = {
    "sampleProcessingProtocol",
    "dataProcessingProtocol",
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "pride_full_scan.log"),
    ],
)
log = logging.getLogger("pride_full_scan")
_last_req = 0.0

PXD_RE = re.compile(r"\b(PXD\d+)\b", re.I)
GEO_RE = re.compile(r"\b(GSE\d+)\b", re.I)
PRJ_RE = re.compile(r"\b(PRJ[NED][A-Z]?\d+)\b", re.I)
PMID_RE = re.compile(r"^\d+$")


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def throttle() -> None:
    global _last_req
    elapsed = time.monotonic() - _last_req
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_req = time.monotonic()


def pride_get(url: str, params: dict | None = None) -> Any:
    last_err: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        throttle()
        try:
            r = http_get(url, params=params, timeout=TIMEOUT)
            if r.status_code == 429 or r.status_code >= 500:
                wait = min(2 ** attempt, 30)
                log.warning("HTTP %s %s (attempt %s), sleep %ss", r.status_code, url, attempt, wait)
                time.sleep(wait)
                continue
            if r.status_code == 404:
                return None
            r.raise_for_status()
            # list/files return JSON; some detail endpoints claim text/plain but are JSON
            ctype = r.headers.get("content-type", "")
            if "json" in ctype or url.rstrip("/").endswith("/projects") or "/files" in url:
                return r.json(), dict(r.headers)
            try:
                return r.json(), dict(r.headers)
            except Exception:
                return r.text, dict(r.headers)
        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 30)
            log.warning("Request failed attempt %s %s: %s; sleep %ss", attempt, url, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"PRIDE request failed after {MAX_RETRIES}: {last_err}")


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
        log.info("Created staging DB %s from %s", STAGING_DB, SCHEMA)
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS src_pride_project (
            accession              TEXT PRIMARY KEY,  -- PXD...
            title                  TEXT,
            description            TEXT,
            organisms              TEXT,   -- JSON array of names
            organism_parts         TEXT,   -- JSON array
            diseases               TEXT,   -- JSON array
            experiment_types       TEXT,   -- JSON array
            instruments            TEXT,   -- JSON array
            softwares              TEXT,   -- JSON array
            keywords               TEXT,   -- JSON array
            project_tags           TEXT,   -- JSON array
            doi                    TEXT,
            license                TEXT,
            submission_date        TEXT,
            publication_date       TEXT,
            submission_type        TEXT,
            related_omics          TEXT,   -- JSON array of related accessions/links
            references_json        TEXT,   -- JSON array of {pubmedID,doi,referenceLine}
            files                  TEXT,   -- JSON compact file list
            file_count             INTEGER,
            raw_json               TEXT,   -- compact project snapshot (no protocols)
            fetched_at             TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_src_pride_fetched_at ON src_pride_project(fetched_at);
        CREATE INDEX IF NOT EXISTS idx_src_pride_pubdate ON src_pride_project(publication_date);
        CREATE TABLE IF NOT EXISTS pride_merge_state (
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
            "page": 0,
            "page_size": PAGE_SIZE,
            "total": None,
            "scanned": 0,
            "skipped_existing": 0,
            "ingested_new": 0,
            "src_upserted": 0,
            "artifacts_created": 0,
            "publications_linked": 0,
            "errors": 0,
            "files_fetched_projects": 0,
        }
    cur = json.loads(row[0])
    cur.setdefault("page_size", PAGE_SIZE)
    return cur


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


def cv_names(items: Any) -> list[str]:
    out: list[str] = []
    if not items:
        return out
    for it in items:
        if isinstance(it, dict):
            name = it.get("name") or it.get("value")
            if name:
                out.append(str(name))
        elif isinstance(it, str):
            out.append(it)
    return out


def join_names(items: Any, sep: str = "; ") -> str | None:
    names = cv_names(items)
    return sep.join(names) if names else None


def parse_related_omics(links: list | None) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for raw in links or []:
        s = str(raw)
        for m in PXD_RE.finditer(s):
            acc = m.group(1).upper()
            if acc not in seen:
                seen.add(acc)
                out.append(acc)
        for m in GEO_RE.finditer(s):
            g = m.group(1).upper()
            if g not in seen:
                seen.add(g)
                out.append(g)
        for m in PRJ_RE.finditer(s):
            p = m.group(1).upper()
            if p not in seen:
                seen.add(p)
                out.append(p)
    return out


def guess_file_kind_role(category: str | None, filename: str) -> tuple[str, str]:
    cat = (category or "").upper()
    name = (filename or "").lower()
    if cat in ("RAW",) or name.endswith((".raw", ".wiff", ".d", ".mzml", ".mzxml")):
        return "raw", "raw_ms"
    if cat in ("PEAK", "SPECTRUM") or name.endswith((".mgf", ".pkl", ".ms2")):
        return "processed", "peak_list"
    if cat in ("SEARCH", "RESULT", "FASTA") or name.endswith(
        (".dat", ".mzid", ".mzTab", ".mztab", ".pep.xml", ".prot.xml")
    ):
        return "processed", "search_result"
    if cat in ("QUANT",) or "ibaq" in name or "lfq" in name:
        return "processed", "quant"
    if cat in ("OTHER", "EXPERIMENTAL_DESIGN") or name.endswith((".txt", ".tsv", ".csv", ".xlsx")):
        return "supplementary", "other"
    return "supplementary", (category or "other").lower()


def format_from_name(name: str) -> str | None:
    if "." not in name:
        return None
    return name.rsplit(".", 1)[-1].lower()


def compact_project(p: dict) -> dict:
    out = {k: v for k, v in p.items() if k not in RAW_DROP_KEYS}
    # trim very large nested blobs defensively
    for k in ("additionalAttributes", "sampleAttributes", "identifiedPTMStrings"):
        if isinstance(out.get(k), list) and len(out[k]) > 50:
            out[k] = out[k][:50]
    return out


def fetch_project_files(accession: str) -> list[dict]:
    """Paginated file metadata; capped and compact."""
    files: list[dict] = []
    page = 0
    total: int | None = None
    while len(files) < MAX_FILES_PER_PROJECT:
        data, headers = pride_get(
            f"{API}/projects/{accession}/files",
            params={"pageSize": FILE_PAGE_SIZE, "page": page},
        )
        if data is None:
            break
        if not isinstance(data, list):
            log.warning("Unexpected files payload for %s: %s", accession, type(data))
            break
        if total is None:
            try:
                total = int(headers.get("total_records") or 0)
            except (TypeError, ValueError):
                total = None
        if not data:
            break
        for f in data:
            if len(files) >= MAX_FILES_PER_PROJECT:
                break
            if not isinstance(f, dict):
                continue
            cat = f.get("fileCategory") or {}
            locs = f.get("publicFileLocations") or []
            ftp = None
            aspera = None
            for L in locs:
                if not isinstance(L, dict):
                    continue
                val = L.get("value") or ""
                name = (L.get("name") or "").lower()
                if "ftp" in name or val.startswith("ftp://"):
                    ftp = val
                elif "aspera" in name or "fasp" in val:
                    aspera = val
            files.append(
                {
                    "name": f.get("fileName"),
                    "cat": cat.get("value") or cat.get("name"),
                    "sz": f.get("fileSizeBytes"),
                    "ftp": ftp,
                    "aspera": aspera,
                    "checksum": f.get("checksum") or None,
                }
            )
        page += 1
        if total is not None and page * FILE_PAGE_SIZE >= total:
            break
        if len(data) < FILE_PAGE_SIZE:
            break
    return files


def upsert_publication(conn: sqlite3.Connection, ref: dict, fetched_at: str) -> str | None:
    pmid = str(ref.get("pubmedID") or "").strip()
    doi = (ref.get("doi") or "").strip() or None
    line = ref.get("referenceLine") or ""
    if pmid and PMID_RE.match(pmid):
        pk = pmid
    elif doi:
        pk = f"DOI:{doi}"
    else:
        return None
    conn.execute(
        """
        INSERT INTO publications(pmid, doi, title, authors, journal, year, abstract, url, fetched_at)
        VALUES (?, ?, ?, NULL, NULL, NULL, NULL, ?, ?)
        ON CONFLICT(pmid) DO UPDATE SET
            doi=COALESCE(excluded.doi, publications.doi),
            title=COALESCE(excluded.title, publications.title),
            fetched_at=excluded.fetched_at
        """,
        (
            pk,
            doi,
            line[:500] if line else None,
            f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if PMID_RE.match(pmid or "") else (
                f"https://doi.org/{doi}" if doi else None
            ),
            fetched_at,
        ),
    )
    return pk


def staging_ingest_project(conn: sqlite3.Connection, project: dict, cur: dict) -> None:
    accession = (project.get("accession") or "").strip().upper()
    if not accession.startswith("PXD"):
        cur["errors"] = cur.get("errors", 0) + 1
        log.warning("skip non-PXD accession: %r", project.get("accession"))
        return

    exists = conn.execute(
        "SELECT 1 FROM src_pride_project WHERE accession=? LIMIT 1", (accession,)
    ).fetchone()
    if exists:
        cur["skipped_existing"] = cur.get("skipped_existing", 0) + 1
        return

    fetched_at = now_iso()
    title = project.get("title") or accession
    description = project.get("projectDescription")
    if isinstance(description, str) and len(description) > 8000:
        description = description[:8000] + "…"

    organisms = cv_names(project.get("organisms"))
    organism_parts = cv_names(project.get("organismParts"))
    diseases = cv_names(project.get("diseases"))
    experiment_types = cv_names(project.get("experimentTypes"))
    instruments = cv_names(project.get("instruments"))
    softwares = cv_names(project.get("softwares"))
    keywords = [str(x) for x in (project.get("keywords") or []) if x]
    project_tags = [str(x) for x in (project.get("projectTags") or []) if x]
    related = parse_related_omics(project.get("otherOmicsLinks"))
    related = [x for x in related if x != accession]

    refs_raw = project.get("references") or []
    refs_compact = []
    for ref in refs_raw:
        if not isinstance(ref, dict):
            continue
        refs_compact.append(
            {
                "pubmedID": ref.get("pubmedID"),
                "doi": ref.get("doi"),
                "referenceLine": (ref.get("referenceLine") or "")[:400],
            }
        )

    files_json: list[dict] = []
    if FETCH_FILES:
        try:
            files_json = fetch_project_files(accession)
            cur["files_fetched_projects"] = cur.get("files_fetched_projects", 0) + 1
        except Exception as e:
            log.warning("files fetch failed %s: %s", accession, e)
            cur["errors"] = cur.get("errors", 0) + 1

    raw_compact = compact_project(project)
    # drop submitter emails from stored snapshot
    if isinstance(raw_compact.get("submitters"), list):
        cleaned = []
        for s in raw_compact["submitters"][:20]:
            if isinstance(s, dict):
                s = {k: v for k, v in s.items() if k != "email"}
            cleaned.append(s)
        raw_compact["submitters"] = cleaned
    if isinstance(raw_compact.get("labPIs"), list):
        cleaned = []
        for s in raw_compact["labPIs"][:20]:
            if isinstance(s, dict):
                s = {k: v for k, v in s.items() if k != "email"}
            cleaned.append(s)
        raw_compact["labPIs"] = cleaned

    conn.execute(
        """
        INSERT INTO src_pride_project(
            accession, title, description, organisms, organism_parts, diseases,
            experiment_types, instruments, softwares, keywords, project_tags,
            doi, license, submission_date, publication_date, submission_type,
            related_omics, references_json, files, file_count, raw_json, fetched_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(accession) DO UPDATE SET
            title=excluded.title, description=excluded.description,
            organisms=excluded.organisms, organism_parts=excluded.organism_parts,
            diseases=excluded.diseases, experiment_types=excluded.experiment_types,
            instruments=excluded.instruments, softwares=excluded.softwares,
            keywords=excluded.keywords, project_tags=excluded.project_tags,
            doi=excluded.doi, license=excluded.license,
            submission_date=excluded.submission_date,
            publication_date=excluded.publication_date,
            submission_type=excluded.submission_type,
            related_omics=excluded.related_omics,
            references_json=excluded.references_json,
            files=excluded.files, file_count=excluded.file_count,
            raw_json=excluded.raw_json, fetched_at=excluded.fetched_at
        """,
        (
            accession,
            title,
            description,
            json.dumps(organisms, ensure_ascii=False),
            json.dumps(organism_parts, ensure_ascii=False),
            json.dumps(diseases, ensure_ascii=False),
            json.dumps(experiment_types, ensure_ascii=False),
            json.dumps(instruments, ensure_ascii=False),
            json.dumps(softwares, ensure_ascii=False),
            json.dumps(keywords, ensure_ascii=False),
            json.dumps(project_tags, ensure_ascii=False),
            project.get("doi"),
            project.get("license"),
            project.get("submissionDate"),
            project.get("publicationDate"),
            project.get("submissionType"),
            json.dumps(related, ensure_ascii=False),
            json.dumps(refs_compact, ensure_ascii=False),
            json.dumps(files_json, ensure_ascii=False),
            len(files_json),
            json.dumps(raw_compact, ensure_ascii=False, default=str),
            fetched_at,
        ),
    )
    cur["src_upserted"] = cur.get("src_upserted", 0) + 1

    canonical_id = f"pride:{accession}"
    ch = content_hash(
        {
            "acc": accession,
            "title": title,
            "org": organisms,
            "assay": experiment_types,
            "nfiles": len(files_json),
        }
    )
    existing_ds = conn.execute(
        "SELECT first_seen_at FROM datasets WHERE canonical_id=?", (canonical_id,)
    ).fetchone()
    first_seen = existing_ds[0] if existing_ds else fetched_at
    organism = join_names(project.get("organisms"))
    tissue = join_names(project.get("organismParts"))
    disease_tags = json.dumps(diseases, ensure_ascii=False) if diseases else None
    assay = join_names(project.get("experimentTypes")) or "proteomics"
    platforms = json.dumps(instruments, ensure_ascii=False) if instruments else None

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
            license=COALESCE(excluded.license, datasets.license),
            content_hash=excluded.content_hash, updated_at=excluded.updated_at,
            probe_status=excluded.probe_status, last_checked=excluded.last_checked
        """,
        (
            canonical_id,
            title,
            title_norm(title),
            organism,
            tissue,
            disease_tags,
            assay,
            platforms,
            None,  # PRIDE list has no reliable sample count
            "public",
            project.get("license"),
            "metadata_only",
            ch,
            first_seen,
            fetched_at,
            "staging",
            fetched_at,
        ),
    )
    if not existing_ds:
        cur["ingested_new"] = cur.get("ingested_new", 0) + 1

    conn.execute(
        """
        INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
        VALUES (?, 'pride', ?, 'project', ?)
        ON CONFLICT(source, accession) DO UPDATE SET
            dataset_id=excluded.dataset_id, role=excluded.role, external_url=excluded.external_url
        """,
        (
            canonical_id,
            accession,
            f"https://www.ebi.ac.uk/pride/archive/projects/{accession}",
        ),
    )

    # Related PXD kept in src_pride_project.related_omics only — do NOT insert as
    # dataset_members(source='pride') or UNIQUE(source, accession) would collide
    # with the real project's primary member row.
    for rel in related:
        if rel.startswith("PXD"):
            continue
        elif rel.startswith("GSE"):
            conn.execute(
                """
                INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
                VALUES (?, 'geo', ?, 'series_link', ?)
                ON CONFLICT(source, accession) DO NOTHING
                """,
                (
                    canonical_id,
                    rel,
                    f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={rel}",
                ),
            )
        elif rel.startswith("PRJ"):
            conn.execute(
                """
                INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
                VALUES (?, 'sra', ?, 'bioproject_link', ?)
                ON CONFLICT(source, accession) DO NOTHING
                """,
                (
                    canonical_id,
                    rel,
                    f"https://www.ncbi.nlm.nih.gov/bioproject/{rel}",
                ),
            )

    for ref in refs_compact:
        pk = upsert_publication(conn, ref, fetched_at)
        if not pk:
            continue
        conn.execute(
            """
            INSERT INTO dataset_publications(dataset_id, pmid, relation, source)
            VALUES (?, ?, 'primary', 'pride')
            ON CONFLICT(dataset_id, pmid) DO NOTHING
            """,
            (canonical_id, pk),
        )
        cur["publications_linked"] = cur.get("publications_linked", 0) + 1

    conn.execute("DELETE FROM dataset_artifacts WHERE dataset_id=?", (canonical_id,))
    for fi in files_json:
        name = fi.get("name") or "file"
        kind, role = guess_file_kind_role(fi.get("cat"), name)
        uri = fi.get("ftp") or fi.get("aspera")
        if not uri:
            uri = (
                f"https://www.ebi.ac.uk/pride/archive/projects/{accession}/files"
                f"?filename={quote(name)}"
            )
        artifact_id = f"pride:{accession}:{hashlib.sha1(name.encode()).hexdigest()[:16]}"
        ach = content_hash({"uri": uri, "size": fi.get("sz"), "name": name})
        size = fi.get("sz")
        try:
            size_i = int(size) if size is not None else None
        except (TypeError, ValueError):
            size_i = None
        conn.execute(
            """
            INSERT INTO dataset_artifacts(
                artifact_id, dataset_id, sample_id, kind, role, title, filename,
                format, uri, size_bytes, access, license, source, derived_from,
                probe_status, last_checked, content_hash
            ) VALUES (?, ?, NULL, ?, ?, ?, ?, ?, ?, ?, ?, NULL, 'pride', NULL, ?, ?, ?)
            ON CONFLICT(artifact_id) DO UPDATE SET
                dataset_id=excluded.dataset_id, uri=excluded.uri, size_bytes=excluded.size_bytes,
                kind=excluded.kind, role=excluded.role, last_checked=excluded.last_checked,
                content_hash=excluded.content_hash
            """,
            (
                artifact_id,
                canonical_id,
                kind,
                role,
                name,
                name,
                format_from_name(name),
                uri,
                size_i,
                "public",
                "listed",
                fetched_at,
                ach,
            ),
        )
        cur["artifacts_created"] = cur.get("artifacts_created", 0) + 1

    conn.commit()


def fetch_projects_page(page: int, page_size: int) -> tuple[list[dict], int | None]:
    data, headers = pride_get(
        f"{API}/projects",
        params={"pageSize": page_size, "page": page},
    )
    total = None
    try:
        total = int(headers.get("total_records") or 0)
    except (TypeError, ValueError):
        pass
    if data is None:
        return [], total
    if not isinstance(data, list):
        raise RuntimeError(f"Unexpected /projects payload type: {type(data)}")
    return data, total


def main() -> int:
    init_staging()
    conn = connect(STAGING_DB)
    cur = load_cursor(conn)
    page = int(cur.get("page") or 0)
    page_size = int(cur.get("page_size") or PAGE_SIZE)
    log.info(
        "PRIDE staging scan start page=%s page_size=%s db=%s fetch_files=%s",
        page,
        page_size,
        STAGING_DB,
        FETCH_FILES,
    )
    save_cursor(conn, cur)

    try:
        while True:
            try:
                projects, total = fetch_projects_page(page, page_size)
            except Exception as e:
                cur["errors"] = cur.get("errors", 0) + 1
                save_cursor(conn, cur, error=str(e)[:500])
                log.exception("page fetch failed page=%s", page)
                time.sleep(5)
                continue

            if total is not None:
                cur["total"] = total
            if not projects:
                cur["completed_at"] = now_iso()
                log.info(
                    "PRIDE scan complete at page=%s scanned=%s total=%s ingested_new=%s errors=%s",
                    page,
                    cur.get("scanned"),
                    cur.get("total"),
                    cur.get("ingested_new"),
                    cur.get("errors"),
                )
                save_cursor(conn, cur)
                break

            for project in projects:
                cur["scanned"] = cur.get("scanned", 0) + 1
                try:
                    staging_ingest_project(conn, project, cur)
                except Exception as e:
                    cur["errors"] = cur.get("errors", 0) + 1
                    acc = (project or {}).get("accession")
                    log.warning("ingest failed %s: %s", acc, e)
                    conn.rollback()

            page += 1
            cur["page"] = page
            save_cursor(conn, cur)
            if page % 5 == 0 or page <= 2:
                log.info(
                    "progress page=%s/%s scanned=%s new=%s skipped=%s artifacts=%s files_proj=%s errors=%s",
                    page,
                    ((cur.get("total") or 0) + page_size - 1) // page_size if cur.get("total") else "?",
                    cur.get("scanned"),
                    cur.get("ingested_new"),
                    cur.get("skipped_existing"),
                    cur.get("artifacts_created"),
                    cur.get("files_fetched_projects"),
                    cur.get("errors"),
                )
    except KeyboardInterrupt:
        log.info("interrupted; saving cursor page=%s", page)
        cur["page"] = page
        save_cursor(conn, cur, error="interrupted")
        return 130
    finally:
        conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
