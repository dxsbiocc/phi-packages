#!/usr/bin/env python3
"""
MetaboLights full ingest (scan phase).

Writes ONLY to omics_metabolights_staging.db — no writes to omics_v1.db.
Resumable via sync_cursors (source=metabolights_full_scan).
Merge to main is a separate later step (not implemented here).
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
STAGING_DB = _DATA_DIR / "omics_metabolights_staging.db"
SCHEMA = SCHEMA_PATH
API = "https://www.ebi.ac.uk/metabolights/ws"
CURSOR_SOURCE = "metabolights_full_scan"
MIN_INTERVAL = 0.30
MAX_RETRIES = 8
TZ = timezone(timedelta(hours=8))
SAVE_EVERY = 25
LOG_EVERY = 50
BATCH_COMMIT = 10

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "metabolights_full_scan.log"),
    ],
)
log = logging.getLogger("metabolights_full_scan")
_last_req = 0.0

PMID_RE = re.compile(r"^\d+$")


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def throttle() -> None:
    global _last_req
    elapsed = time.monotonic() - _last_req
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_req = time.monotonic()


def mtbls_get(url: str, params: dict | None = None) -> Any:
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
            if r.status_code in (400, 401, 404):
                # 401 = embargoed/private; 400 = bad id; 404 = missing
                return None
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 30)
            log.warning("Request failed attempt %s %s: %s; sleep %ss", attempt, url, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"MetaboLights request failed after {MAX_RETRIES}: {last_err}")


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
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS src_metabolights_study (
            accession              TEXT PRIMARY KEY,  -- MTBLS...
            title                  TEXT,
            organism               TEXT,
            organism_part          TEXT,
            assay                  TEXT,
            technology             TEXT,
            platforms              TEXT,   -- JSON array
            study_status           TEXT,
            sample_count           INTEGER,
            study_size_bytes       INTEGER,
            public_release_date    TEXT,
            submission_date        TEXT,
            descriptors            TEXT,   -- JSON array
            publications           TEXT,   -- JSON array of {pmid,doi,title}
            files                  TEXT,   -- JSON array
            raw_json               TEXT,   -- compact snapshot (no sample rows)
            fetched_at             TEXT
        )
        """
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_src_mtbls_fetched_at "
        "ON src_metabolights_study(fetched_at)"
    )
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS mtbls_merge_state (
            accession TEXT PRIMARY KEY,
            merged_at TEXT
        )
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
            "skipped_private": 0,
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


def ms_to_date(val: Any) -> str | None:
    if val is None:
        return None
    if isinstance(val, str) and re.match(r"^\d{4}-\d{2}-\d{2}", val):
        return val[:10]
    try:
        ms = int(val)
        if ms > 10_000_000_000:  # ms
            ms = ms // 1000
        return datetime.fromtimestamp(ms, tz=timezone.utc).strftime("%Y-%m-%d")
    except (TypeError, ValueError, OSError):
        return str(val)[:32] if val else None


def extract_organisms(org_list: Any) -> tuple[str | None, str | None]:
    if not org_list:
        return None, None
    names: list[str] = []
    parts: list[str] = []
    if isinstance(org_list, list):
        for o in org_list:
            if isinstance(o, dict):
                n = o.get("organismName") or o.get("name")
                p = o.get("organismPart") or o.get("part")
                if n:
                    names.append(str(n))
                if p:
                    parts.append(str(p))
            elif o:
                names.append(str(o))
    elif isinstance(org_list, str):
        names.append(org_list)
    org = "; ".join(dict.fromkeys(names)) if names else None
    part = "; ".join(dict.fromkeys(parts)) if parts else None
    return org, part


def slim_public_study(content: dict) -> dict:
    """Keep useful metadata; drop bulky sample/assay table cells."""
    assays_slim = []
    for a in content.get("assays") or []:
        if not isinstance(a, dict):
            continue
        assays_slim.append(
            {
                "measurement": a.get("measurement"),
                "technology": a.get("technology"),
                "platform": a.get("platform"),
                "fileName": a.get("fileName"),
                "assayNumber": a.get("assayNumber"),
            }
        )
    pubs = []
    for p in content.get("publications") or []:
        if not isinstance(p, dict):
            continue
        pubs.append(
            {
                "title": p.get("title"),
                "doi": p.get("doi"),
                "pubmedId": p.get("pubmedId"),
                "authorList": p.get("authorList"),
            }
        )
    st = content.get("sampleTable") or {}
    sample_n = len(st.get("data") or []) if isinstance(st, dict) else None
    return {
        "studyIdentifier": content.get("studyIdentifier") or content.get("id"),
        "title": content.get("title"),
        "description": (content.get("description") or content.get("studyDescription") or "")[:4000],
        "organism": content.get("organism"),
        "studyStatus": content.get("studyStatus"),
        "studyPublicReleaseDate": content.get("studyPublicReleaseDate"),
        "studySubmissionDate": content.get("studySubmissionDate"),
        "updateDate": content.get("updateDate"),
        "studySize": content.get("studySize"),
        "publicStudy": content.get("publicStudy"),
        "factors": content.get("factors"),
        "descriptors": content.get("descriptors"),
        "publications": pubs,
        "assays": assays_slim,
        "contacts_n": len(content.get("contacts") or []),
        "sample_count": sample_n,
        "protocols_n": len(content.get("protocols") or []),
    }


def fetch_accession_list() -> list[str]:
    data = mtbls_get(f"{API}/studies")
    if not data:
        raise RuntimeError("Failed to fetch MetaboLights study list")
    content = data.get("content") if isinstance(data, dict) else data
    if not isinstance(content, list):
        # some responses nest under studies
        content = (data.get("studies") if isinstance(data, dict) else None) or []
    accs = [str(x).strip() for x in content if x and str(x).startswith("MTBLS")]
    # stable order: numeric suffix
    def sort_key(a: str) -> int:
        m = re.search(r"(\d+)$", a)
        return int(m.group(1)) if m else 0

    accs = sorted(dict.fromkeys(accs), key=sort_key)
    return accs


def staging_ingest(conn: sqlite3.Connection, accession: str, cur: dict) -> None:
    exists = conn.execute(
        "SELECT 1 FROM src_metabolights_study WHERE accession=? LIMIT 1", (accession,)
    ).fetchone()
    if exists:
        cur["skipped_existing"] = cur.get("skipped_existing", 0) + 1
        return

    fetched_at = now_iso()
    public = mtbls_get(f"{API}/studies/public/study/{accession}")
    if not public or not public.get("content"):
        # try ISA endpoint as fallback for public-but-odd studies
        isa = mtbls_get(f"{API}/studies/{accession}")
        if not isa:
            cur["skipped_private"] = cur.get("skipped_private", 0) + 1
            log.warning("skip unavailable %s", accession)
            return
        inv = isa.get("isaInvestigation") or {}
        mt = isa.get("mtblsStudy") or {}
        studies = inv.get("studies") or []
        st0 = studies[0] if studies else {}
        title = inv.get("title") or st0.get("title") or accession
        assays = st0.get("assays") or []
        assay_meas = []
        techs = []
        platforms = []
        for a in assays:
            if not isinstance(a, dict):
                continue
            mt_ = a.get("measurementType") or {}
            tt = a.get("technologyType") or {}
            if isinstance(mt_, dict) and mt_.get("annotationValue"):
                assay_meas.append(mt_["annotationValue"])
            if isinstance(tt, dict) and tt.get("annotationValue"):
                techs.append(tt["annotationValue"])
            if a.get("technologyPlatform"):
                platforms.append(str(a["technologyPlatform"]))
        files_resp = mtbls_get(f"{API}/studies/{accession}/files") or {}
        files_items = files_resp.get("study") or []
        files_json = [
            {
                "file": f.get("file"),
                "type": f.get("type"),
                "directory": f.get("directory"),
                "status": f.get("status"),
                "timestamp": f.get("timestamp"),
            }
            for f in files_items
            if isinstance(f, dict)
        ]
        raw = {
            "source": "isa",
            "title": title,
            "description": (inv.get("description") or "")[:4000],
            "studyStatus": mt.get("studyStatus"),
            "assays": [{"measurement": m} for m in assay_meas],
            "technology": techs,
            "platforms": platforms,
        }
        organism = None
        organism_part = None
        sample_count = len(st0.get("samples") or []) or None
        assay = "; ".join(dict.fromkeys(assay_meas)) or None
        technology = "; ".join(dict.fromkeys(techs)) or None
        status = mt.get("studyStatus")
        pub_date = (inv.get("publicReleaseDate") or "")[:10] or None
        sub_date = (inv.get("submissionDate") or "")[:10] or None
        study_size = None
        descriptors = [
            d.get("annotationValue")
            for d in (st0.get("studyDesignDescriptors") or [])
            if isinstance(d, dict) and d.get("annotationValue")
        ]
        pubs = []
        for p in (st0.get("publications") or inv.get("publications") or []):
            if isinstance(p, dict):
                pubs.append(
                    {
                        "title": p.get("title"),
                        "doi": p.get("doi"),
                        "pubmedId": p.get("pubMedID") or p.get("pubmedId"),
                    }
                )
    else:
        content = public["content"]
        slim = slim_public_study(content)
        title = slim.get("title") or accession
        organism, organism_part = extract_organisms(slim.get("organism"))
        assay_meas = [a.get("measurement") for a in slim.get("assays") or [] if a.get("measurement")]
        techs = [a.get("technology") for a in slim.get("assays") or [] if a.get("technology")]
        platforms = [a.get("platform") for a in slim.get("assays") or [] if a.get("platform")]
        assay = "; ".join(dict.fromkeys(str(x) for x in assay_meas)) or None
        technology = "; ".join(dict.fromkeys(str(x) for x in techs)) or None
        sample_count = slim.get("sample_count")
        status = slim.get("studyStatus")
        pub_date = ms_to_date(slim.get("studyPublicReleaseDate"))
        sub_date = ms_to_date(slim.get("studySubmissionDate"))
        study_size = slim.get("studySize")
        try:
            study_size = int(study_size) if study_size is not None else None
        except (TypeError, ValueError):
            study_size = None
        descriptors = [
            d.get("description")
            for d in (slim.get("descriptors") or [])
            if isinstance(d, dict) and d.get("description")
        ]
        pubs = slim.get("publications") or []
        files_resp = mtbls_get(f"{API}/studies/{accession}/files") or {}
        files_items = files_resp.get("study") or []
        files_json = [
            {
                "file": f.get("file"),
                "type": f.get("type"),
                "directory": f.get("directory"),
                "status": f.get("status"),
                "timestamp": f.get("timestamp"),
            }
            for f in files_items
            if isinstance(f, dict)
        ]
        raw = slim
        raw["source"] = "public_study"

    conn.execute(
        """
        INSERT INTO src_metabolights_study(
            accession, title, organism, organism_part, assay, technology, platforms,
            study_status, sample_count, study_size_bytes, public_release_date,
            submission_date, descriptors, publications, files, raw_json, fetched_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(accession) DO UPDATE SET
            title=excluded.title, organism=excluded.organism,
            organism_part=excluded.organism_part, assay=excluded.assay,
            technology=excluded.technology, platforms=excluded.platforms,
            study_status=excluded.study_status, sample_count=excluded.sample_count,
            study_size_bytes=excluded.study_size_bytes,
            public_release_date=excluded.public_release_date,
            submission_date=excluded.submission_date,
            descriptors=excluded.descriptors, publications=excluded.publications,
            files=excluded.files, raw_json=excluded.raw_json, fetched_at=excluded.fetched_at
        """,
        (
            accession,
            title,
            organism,
            organism_part,
            assay,
            technology,
            json.dumps(platforms) if platforms else None,
            status,
            sample_count,
            study_size,
            pub_date,
            sub_date,
            json.dumps(descriptors) if descriptors else None,
            json.dumps(pubs) if pubs else None,
            json.dumps(files_json),
            json.dumps(raw, ensure_ascii=False),
            fetched_at,
        ),
    )
    cur["src_upserted"] = cur.get("src_upserted", 0) + 1

    canonical_id = f"mtbls:{accession}"
    ch = content_hash({"acc": accession, "title": title, "n": sample_count, "org": organism})
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
            title=excluded.title, title_norm=excluded.title_norm,
            organism=COALESCE(excluded.organism, datasets.organism),
            tissue=COALESCE(excluded.tissue, datasets.tissue),
            assay=COALESCE(excluded.assay, datasets.assay),
            platforms=COALESCE(excluded.platforms, datasets.platforms),
            n_samples=COALESCE(excluded.n_samples, datasets.n_samples),
            content_hash=excluded.content_hash, updated_at=excluded.updated_at,
            probe_status=excluded.probe_status, last_checked=excluded.last_checked
        """,
        (
            canonical_id,
            title,
            title_norm(title),
            organism,
            organism_part,
            None,
            assay or technology,
            json.dumps(platforms) if platforms else None,
            sample_count,
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

    conn.execute(
        """
        INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
        VALUES (?, 'metabolights', ?, 'study', ?)
        ON CONFLICT(source, accession) DO UPDATE SET
            dataset_id=excluded.dataset_id, role=excluded.role, external_url=excluded.external_url
        """,
        (
            canonical_id,
            accession,
            f"https://www.ebi.ac.uk/metabolights/{accession}",
        ),
    )

    # publications (light); tolerate pmid/doi uniqueness clashes
    for p in pubs or []:
        if not isinstance(p, dict):
            continue
        pmid = str(p.get("pubmedId") or "").strip()
        doi = (p.get("doi") or "").strip() or None
        if pmid and PMID_RE.match(pmid):
            pk = pmid
        elif doi:
            pk = f"DOI:{doi}"
        else:
            continue
        # If this DOI already lives under another pmid, reuse that row.
        if doi:
            existing_doi = conn.execute(
                "SELECT pmid FROM publications WHERE doi=? LIMIT 1", (doi,)
            ).fetchone()
            if existing_doi and existing_doi[0] != pk:
                pk = existing_doi[0]
        try:
            conn.execute(
                """
                INSERT INTO publications(pmid, doi, title, authors, journal, year, abstract, url, fetched_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(pmid) DO UPDATE SET
                    title=COALESCE(excluded.title, publications.title),
                    authors=COALESCE(excluded.authors, publications.authors),
                    fetched_at=excluded.fetched_at
                """,
                (
                    pk,
                    doi,
                    p.get("title"),
                    p.get("authorList"),
                    None,
                    None,
                    None,
                    f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if PMID_RE.match(pmid or "") else None,
                    fetched_at,
                ),
            )
        except sqlite3.IntegrityError:
            # DOI unique clash or similar — link to whatever pmid already holds the DOI
            if doi:
                row = conn.execute(
                    "SELECT pmid FROM publications WHERE doi=? LIMIT 1", (doi,)
                ).fetchone()
                if row:
                    pk = row[0]
                else:
                    continue
            else:
                continue
        conn.execute(
            """
            INSERT INTO dataset_publications(dataset_id, pmid, relation, source)
            VALUES (?, ?, 'primary', 'metabolights')
            ON CONFLICT(dataset_id, pmid) DO NOTHING
            """,
            (canonical_id, pk),
        )

    conn.execute("DELETE FROM dataset_artifacts WHERE dataset_id=?", (canonical_id,))
    for fi in files_json:
        name = fi.get("file") or "file"
        if fi.get("directory"):
            continue
        ftype = (fi.get("type") or "").lower()
        if "raw" in ftype:
            kind, role = "raw", "raw"
        elif "maf" in ftype or "metabolite" in ftype:
            kind, role = "processed", "metabolite_assignment"
        elif "metadata" in ftype or ftype.startswith("metadata"):
            kind, role = "supplementary", "isa_tab"
        else:
            kind, role = "supplementary", ftype or "other"
        uri = (
            f"https://www.ebi.ac.uk/metabolights/{accession}/files/{quote(name, safe='')}"
        )
        artifact_id = f"mtbls:{accession}:{hashlib.sha1(name.encode()).hexdigest()[:16]}"
        ach = content_hash({"uri": uri, "name": name, "type": ftype})
        conn.execute(
            """
            INSERT INTO dataset_artifacts(
                artifact_id, dataset_id, sample_id, kind, role, title, filename,
                format, uri, size_bytes, access, license, source, derived_from,
                probe_status, last_checked, content_hash
            ) VALUES (?, ?, NULL, ?, ?, ?, ?, ?, ?, NULL, 'public', NULL, 'metabolights',
                      NULL, 'staging', ?, ?)
            ON CONFLICT(artifact_id) DO UPDATE SET
                kind=excluded.kind, role=excluded.role, uri=excluded.uri,
                last_checked=excluded.last_checked, content_hash=excluded.content_hash
            """,
            (
                artifact_id,
                canonical_id,
                kind,
                role,
                name,
                name,
                name.rsplit(".", 1)[-1].lower() if "." in name else None,
                uri,
                fetched_at,
                ach,
            ),
        )
        cur["artifacts_created"] = cur.get("artifacts_created", 0) + 1


def main() -> int:
    log.info("Starting MetaboLights SCAN-ONLY (staging=%s)", STAGING_DB)
    init_staging()
    st = connect(STAGING_DB, busy_ms=5000)
    cur = load_cursor(st)

    try:
        log.info("Fetching public study accession list...")
        accs = fetch_accession_list()
        total = len(accs)
        cur["total"] = total
        idx = int(cur.get("index") or 0)
        if idx < 0:
            idx = 0
        if idx > total:
            idx = total
        save_cursor(st, cur)
        log.info("Resume index=%s / total=%s cursor=%s", idx, total, json.dumps(cur))

        since_log = 0
        since_save = 0
        pending_commit = 0

        while idx < total:
            accession = accs[idx]
            try:
                staging_ingest(st, accession, cur)
                pending_commit += 1
            except Exception as e:
                cur["errors"] = cur.get("errors", 0) + 1
                log.exception("Failed staging %s: %s", accession, e)
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
                    "progress index=%s/%s scanned=%s new=%s skipped=%s src=%s artifacts=%s private=%s errors=%s",
                    idx,
                    total,
                    cur.get("scanned"),
                    cur.get("ingested_new"),
                    cur.get("skipped_existing"),
                    cur.get("src_upserted"),
                    cur.get("artifacts_created"),
                    cur.get("skipped_private"),
                    cur.get("errors"),
                )

        if pending_commit:
            st.commit()
        cur["index"] = total
        cur["completed_at"] = now_iso()
        save_cursor(st, cur)
        log.info("Finished MetaboLights scan: %s", json.dumps(cur))
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
