#!/usr/bin/env python3
"""
ArrayExpress / BioStudies full ingest (scan phase).

Writes ONLY to omics_ae_staging.db during the scan — no writes to omics_v1.db —
so sra_orphan_scan can hold the main DB without lock contention.

Resumable via sync_cursors in staging (source=arrayexpress_full_scan).
After scan completes (or anytime offline), run arrayexpress_merge_to_main.py
to drain staging into omics_v1.db with L1 dedup.
"""

from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)

import csv
import hashlib
import io
import json
import logging
import re
import sqlite3
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable
from urllib.parse import quote

import requests

from http_util import TIMEOUT, http_get

from paths import CODE_ROOT as ROOT, SCHEMA_PATH, data_dir as _data_dir, logs_dir as _logs_dir
_DATA_DIR = _data_dir(); _DATA_DIR.mkdir(parents=True, exist_ok=True)
_LOGS_DIR = _logs_dir(); _LOGS_DIR.mkdir(parents=True, exist_ok=True)
MAIN_DB = _DATA_DIR / "omics_v1.db"
STAGING_DB = _DATA_DIR / "omics_ae_staging.db"
SCHEMA = SCHEMA_PATH
API = "https://www.ebi.ac.uk/biostudies/api/v1"
CURSOR_SOURCE = "arrayexpress_full_scan"
PAGE_SIZE = 10  # keep small; EBI search hard-fails past ~offset 20000. Year windows keep each query under that limit.
MIN_INTERVAL = 0.22
MAX_RETRIES = 8
TZ = timezone(timedelta(hours=8))
MERGE_EVERY = 0  # scan-only: never live-merge to main (use arrayexpress_merge_to_main.py)
WRITE_MAIN_DURING_SCAN = False
MAX_SDRF_BYTES = 2_000_000
MAX_SAMPLES_PER_STUDY = 300
SDRF_FETCH = False  # first full pass: metadata+files only; SDRF enrich later
LOCK_RETRY_WAIT = 3.0
LOCK_RETRY_MAX = 60   # wait up to ~3 min for orphan to release write lock

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "arrayexpress_full_scan.log"),
    ],
)
log = logging.getLogger("arrayexpress_full_scan")
_last_req = 0.0

GEO_RE = re.compile(r"\b(GSE\d+)\b", re.I)
PRJ_RE = re.compile(r"\b(PRJ[NED][A-Z]?\d+)\b", re.I)
SRA_STUDY_RE = re.compile(r"\b([SED]RP\d+)\b", re.I)
SRA_EXP_RE = re.compile(r"\b([SED]RX\d+)\b", re.I)
AE_ACC_RE = re.compile(r"^E-[A-Z]+-\d+$", re.I)


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def throttle() -> None:
    global _last_req
    elapsed = time.monotonic() - _last_req
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_req = time.monotonic()


def ebi_get(url: str, params: dict | None = None, raw: bool = False) -> Any:
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
            if raw:
                return r.content
            ctype = r.headers.get("content-type", "")
            if "json" in ctype or "/search" in url or "/studies/" in url:
                return r.json()
            return r.text
        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 30)
            log.warning("Request failed attempt %s %s: %s; sleep %ss", attempt, url, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"EBI request failed after {MAX_RETRIES}: {last_err}")


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
        sql = SCHEMA.read_text()
        conn.executescript(sql)
        conn.commit()
        log.info("Created staging DB %s", STAGING_DB)
    # merge watermark table
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS ae_merge_state (
            accession TEXT PRIMARY KEY,
            merged_at TEXT
        )
        """
    )
    conn.commit()
    conn.close()


def with_lock_retry(fn: Callable[[], Any], what: str) -> Any:
    last: Exception | None = None
    for attempt in range(1, LOCK_RETRY_MAX + 1):
        try:
            return fn()
        except sqlite3.OperationalError as e:
            if "locked" not in str(e).lower():
                raise
            last = e
            if attempt % 6 == 1:
                log.warning("DB locked during %s (attempt %s); retrying...", what, attempt)
            time.sleep(LOCK_RETRY_WAIT)
    raise RuntimeError(f"Gave up on {what}: {last}")


def load_cursor(conn: sqlite3.Connection) -> dict:
    row = conn.execute(
        "SELECT cursor_value FROM sync_cursors WHERE source=?", (CURSOR_SOURCE,)
    ).fetchone()
    if not row or not row[0]:
        return {
            "page": 1,
            "total": None,
            "scanned": 0,
            "skipped_existing": 0,
            "ingested_new": 0,
            "merged": 0,
            "samples_created": 0,
            "artifacts_created": 0,
            "src_upserted": 0,
            "merged_to_main": 0,
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
        (CURSOR_SOURCE, json.dumps(cur), "page", now_iso(), error),
    )
    conn.commit()


def mirror_cursor_to_main(cur: dict, error: str | None = None) -> None:
    """Best-effort cursor mirror; disabled during scan-only mode."""
    if not WRITE_MAIN_DURING_SCAN:
        return
    for attempt in range(1, 4):
        try:
            conn = connect(MAIN_DB, busy_ms=15000)
            try:
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
                return
            finally:
                conn.close()
        except sqlite3.OperationalError as e:
            if "locked" in str(e).lower() and attempt < 3:
                time.sleep(1.5)
                continue
            log.warning("Could not mirror cursor to main: %s", e)
            return
        except Exception as e:
            log.warning("Could not mirror cursor to main: %s", e)
            return


def attr_map(attrs: list | None) -> dict[str, str]:
    out: dict[str, str] = {}
    for a in attrs or []:
        if isinstance(a, dict) and a.get("name") is not None and a.get("name") not in out:
            out[str(a["name"])] = "" if a.get("value") is None else str(a["value"])
    return out


def walk_collect(node: Any, attrs_acc: dict[str, list[str]], links_acc: list[dict]) -> None:
    if isinstance(node, list):
        for x in node:
            walk_collect(x, attrs_acc, links_acc)
        return
    if not isinstance(node, dict):
        return
    for a in node.get("attributes") or []:
        if isinstance(a, dict) and a.get("name") is not None:
            attrs_acc.setdefault(str(a["name"]), []).append(
                "" if a.get("value") is None else str(a["value"])
            )
    for L in node.get("links") or []:
        if isinstance(L, list):
            for x in L:
                if isinstance(x, dict):
                    links_acc.append(x)
        elif isinstance(L, dict):
            links_acc.append(L)
    for sub in node.get("subsections") or []:
        walk_collect(sub, attrs_acc, links_acc)


def link_type(link: dict) -> str:
    for a in link.get("attributes") or []:
        if isinstance(a, dict) and (a.get("name") or "").lower() == "type":
            return str(a.get("value") or "")
    return ""


def extract_related(links: list[dict], text_blobs: list[str]) -> tuple[list[str], list[str], list[str]]:
    geos: list[str] = []
    sras: list[str] = []
    prjs: list[str] = []
    seen_g, seen_s, seen_p = set(), set(), set()

    def add_geo(x: str) -> None:
        x = x.upper()
        if x not in seen_g:
            seen_g.add(x)
            geos.append(x)

    def add_sra(x: str) -> None:
        x = x.upper()
        if x not in seen_s:
            seen_s.add(x)
            sras.append(x)

    def add_prj(x: str) -> None:
        x = x.upper()
        if x not in seen_p:
            seen_p.add(x)
            prjs.append(x)

    for L in links:
        url = str(L.get("url") or "")
        typ = link_type(L).lower()
        if typ == "geo" or GEO_RE.fullmatch(url):
            m = GEO_RE.search(url)
            if m:
                add_geo(m.group(1))
        if typ in ("ena", "sra", "bioproject", "ena study", "ena experiment"):
            for m in PRJ_RE.finditer(url):
                add_prj(m.group(1))
            for m in SRA_STUDY_RE.finditer(url):
                add_sra(m.group(1))
            for m in SRA_EXP_RE.finditer(url):
                add_sra(m.group(1))
            if re.fullmatch(r"[SED]RP\d+", url, re.I):
                add_sra(url.upper())
            if re.fullmatch(r"PRJ[NED][A-Z]?\d+", url, re.I):
                add_prj(url.upper())
        for m in GEO_RE.finditer(url):
            add_geo(m.group(1))
        for m in PRJ_RE.finditer(url):
            add_prj(m.group(1))
        for m in SRA_STUDY_RE.finditer(url):
            add_sra(m.group(1))

    blob = " ".join(text_blobs)
    for m in GEO_RE.finditer(blob):
        add_geo(m.group(1))
    for m in PRJ_RE.finditer(blob):
        add_prj(m.group(1))
    for m in SRA_STUDY_RE.finditer(blob):
        add_sra(m.group(1))
    return geos, sras, prjs


def first_attr(attrs: dict[str, list[str]], *names: str) -> str | None:
    for n in names:
        vals = attrs.get(n)
        if vals:
            return vals[0]
    lower = {k.lower(): v for k, v in attrs.items()}
    for n in names:
        vals = lower.get(n.lower())
        if vals:
            return vals[0]
    return None


def guess_artifact_kind_role(name: str, typ: str, desc: str, section: str) -> tuple[str, str]:
    blob = f"{name} {typ} {desc} {section}".lower()
    if "sdrf" in blob or "idf" in blob or "mage-tab" in blob:
        return "supplementary", "mage_tab"
    if any(x in blob for x in ("processed", "count", "tpm", "fpkm", "matrix", "normalized")):
        role = "matrix"
        if "tpm" in blob:
            role = "tpm"
        elif "fpkm" in blob:
            role = "fpkm"
        elif "count" in blob:
            role = "counts"
        return "processed", role
    if any(x in blob for x in ("raw", "fastq", "bam", "cel")):
        return "raw", "raw"
    return "supplementary", "other"


def format_from_name(name: str) -> str | None:
    if "." not in name:
        return None
    return name.rsplit(".", 1)[-1].lower()


def parse_sdrf_samples(content: bytes, dataset_id: str, accession: str) -> list[dict]:
    try:
        text = content.decode("utf-8", errors="replace")
    except Exception:
        return []
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        return []
    reader = csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t")
    if not reader.fieldnames:
        return []

    def col(*names: str) -> str | None:
        lower = {c.lower(): c for c in reader.fieldnames or []}
        for n in names:
            if n.lower() in lower:
                return lower[n.lower()]
        for n in names:
            for k, orig in lower.items():
                if n.lower() in k:
                    return orig
        return None

    src_col = col("Source Name", "Sample Name", "Assay Name")
    org_col = col("Characteristics[Organism]", "Characteristics [Organism]", "Organism")
    tissue_col = col(
        "Characteristics[Organism part]",
        "Characteristics[organism part]",
        "Characteristics[Tissue]",
        "Characteristics[tissue]",
    )
    disease_col = col("Characteristics[Disease]", "Characteristics[disease]", "Characteristics[Disease State]")
    sex_col = col("Characteristics[Sex]", "Characteristics[sex]", "Characteristics[Gender]")
    age_col = col("Characteristics[Age]", "Characteristics[age]")
    factor_cols = [c for c in (reader.fieldnames or []) if c.lower().startswith("factor value")]

    samples: list[dict] = []
    seen: set[str] = set()
    for row in reader:
        if len(samples) >= MAX_SAMPLES_PER_STUDY:
            break
        name = (row.get(src_col) or "").strip() if src_col else ""
        if not name:
            for v in row.values():
                if v and str(v).strip():
                    name = str(v).strip()
                    break
        if not name or name in seen:
            continue
        seen.add(name)
        attrs = {k: v for k, v in row.items() if v}
        condition = None
        if factor_cols:
            condition = "; ".join(f"{c}={row.get(c)}" for c in factor_cols if row.get(c))
        sample_id = f"ae:{accession}:{hashlib.sha1(name.encode()).hexdigest()[:12]}"
        samples.append(
            {
                "sample_id": sample_id,
                "dataset_id": dataset_id,
                "title": name,
                "organism": (row.get(org_col) if org_col else None) or None,
                "tissue": (row.get(tissue_col) if tissue_col else None) or None,
                "disease": (row.get(disease_col) if disease_col else None) or None,
                "condition": condition,
                "sex": (row.get(sex_col) if sex_col else None) or None,
                "age": (row.get(age_col) if age_col else None) or None,
                "attributes": attrs,
            }
        )
    return samples


def staging_ingest_study(conn: sqlite3.Connection, accession: str, hit: dict, cur: dict) -> None:
    """Write study into staging with provisional canonical_id=ae:<acc>; L1 merge happens on flush to main."""
    exists = conn.execute(
        "SELECT 1 FROM src_arrayexpress_entry WHERE accession=? LIMIT 1", (accession,)
    ).fetchone()
    if exists:
        cur["skipped_existing"] = cur.get("skipped_existing", 0) + 1
        return
    fetched_at = now_iso()
    study = ebi_get(f"{API}/studies/{accession}")
    if not study:
        cur["errors"] = cur.get("errors", 0) + 1
        log.warning("study not found %s", accession)
        return

    attrs_acc: dict[str, list[str]] = {}
    links_acc: list[dict] = []
    walk_collect(study.get("section") or {}, attrs_acc, links_acc)
    top = attr_map(study.get("attributes"))
    title = top.get("Title") or first_attr(attrs_acc, "Title") or hit.get("title") or accession
    organism = first_attr(attrs_acc, "Organism")
    assay = first_attr(attrs_acc, "Study type", "Technology", "Assay by Molecule")
    sample_count_s = first_attr(attrs_acc, "Sample count")
    try:
        sample_count = int(sample_count_s) if sample_count_s else None
    except ValueError:
        sample_count = None

    geos, sras, prjs = extract_related(links_acc, [title or "", json.dumps(attrs_acc)[:5000]])
    m_geod = re.match(r"E-GEOD-(\d+)$", accession, re.I)
    if m_geod:
        gse = f"GSE{m_geod.group(1)}"
        if gse not in geos:
            geos.append(gse)

    related_geo = ",".join(geos) if geos else None
    related_sra = ",".join(sras + prjs) if (sras or prjs) else None

    info = ebi_get(f"{API}/studies/{accession}/info") or {}
    http_base = (info.get("httpLink") or "").rstrip("/")
    ftp_base = (info.get("ftpLink") or "").rstrip("/")

    files_items: list[dict] = []
    offset = 0
    limit = 100
    while True:
        fj = ebi_get(f"{API}/studies/{accession}/files", params={"limit": limit, "offset": offset})
        if not fj:
            break
        items = fj.get("items") or []
        files_items.extend(items)
        pag = fj.get("pagination") or {}
        total_f = int(pag.get("total") or 0)
        offset += len(items)
        if not items or offset >= total_f or offset > 5000:
            break

    files_json = []
    for it in files_items:
        name = it.get("Name") or it.get("path") or ""
        path = it.get("path") or name
        size = it.get("Size")
        try:
            size_i = int(size) if size is not None else None
        except (TypeError, ValueError):
            size_i = None
        files_json.append(
            {
                "name": name,
                "path": path,
                "size": size_i,
                "type": it.get("Type") or it.get("file_type"),
                "description": it.get("Description"),
                "section": it.get("Section"),
                "format": it.get("Format"),
            }
        )

    raw_compact = {
        "accno": study.get("accno"),
        "attributes": study.get("attributes"),
        "hit": {k: hit.get(k) for k in ("accession", "title", "release_date", "files", "links", "isPublic")},
        "links": links_acc[:50],
        "related_geo": geos,
        "related_sra": sras,
        "related_prj": prjs,
        "key_attrs": {k: v[:5] for k, v in list(attrs_acc.items())[:40]},
        "http_base": http_base,
        "ftp_base": ftp_base,
    }

    conn.execute(
        """
        INSERT INTO src_arrayexpress_entry(
            accession, title, organism, assay, sample_count,
            related_geo, related_sra, files, raw_json, fetched_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(accession) DO UPDATE SET
            title=excluded.title, organism=excluded.organism, assay=excluded.assay,
            sample_count=excluded.sample_count, related_geo=excluded.related_geo,
            related_sra=excluded.related_sra, files=excluded.files,
            raw_json=NULL, fetched_at=excluded.fetched_at
        """,
        (
            accession, title, organism, assay, sample_count, related_geo, related_sra,
            json.dumps(files_json), None, fetched_at,  # raw_json omitted
        ),
    )
    cur["src_upserted"] = cur.get("src_upserted", 0) + 1

    # provisional dataset in staging
    canonical_id = f"ae:{accession}"
    ch = content_hash({"acc": accession, "title": title, "n": sample_count})
    existing = conn.execute(
        "SELECT first_seen_at FROM datasets WHERE canonical_id=?", (canonical_id,)
    ).fetchone()
    first_seen = existing[0] if existing else fetched_at
    access = "public" if hit.get("isPublic", True) else "unknown"
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
            n_samples=COALESCE(excluded.n_samples, datasets.n_samples),
            content_hash=excluded.content_hash, updated_at=excluded.updated_at,
            probe_status=excluded.probe_status, last_checked=excluded.last_checked
        """,
        (
            canonical_id, title, title_norm(title), organism, None, None, assay, None,
            sample_count, access, None, "metadata_only", ch, first_seen, fetched_at,
            "staging", fetched_at,
        ),
    )
    if not existing:
        cur["ingested_new"] = cur.get("ingested_new", 0) + 1

    conn.execute(
        """
        INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
        VALUES (?, 'arrayexpress', ?, 'entry', ?)
        ON CONFLICT(source, accession) DO UPDATE SET
            dataset_id=excluded.dataset_id, role=excluded.role, external_url=excluded.external_url
        """,
        (canonical_id, accession, f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}"),
    )
    # store related accessions as members in staging for merge logic
    for g in geos:
        conn.execute(
            """
            INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
            VALUES (?, 'geo', ?, 'series_link', ?)
            ON CONFLICT(source, accession) DO NOTHING
            """,
            (canonical_id, g, f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}"),
        )
    for p in prjs:
        conn.execute(
            """
            INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
            VALUES (?, 'sra', ?, 'bioproject_link', ?)
            ON CONFLICT(source, accession) DO NOTHING
            """,
            (canonical_id, p, f"https://www.ncbi.nlm.nih.gov/bioproject/{p}"),
        )
    for s in sras:
        role = "study_link" if re.match(r"[SED]RP\d+$", s, re.I) else "experiment_link"
        conn.execute(
            """
            INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
            VALUES (?, 'sra', ?, ?, ?)
            ON CONFLICT(source, accession) DO NOTHING
            """,
            (canonical_id, s, role, f"https://www.ebi.ac.uk/ena/browser/view/{s}"),
        )

    # clear previous staging artifacts/samples for this provisional id then reinsert
    conn.execute("DELETE FROM dataset_artifacts WHERE dataset_id=?", (canonical_id,))
    conn.execute("DELETE FROM samples WHERE dataset_id=?", (canonical_id,))

    for fi in files_json:
        name = fi.get("name") or fi.get("path") or "file"
        path = fi.get("path") or name
        kind, role = guess_artifact_kind_role(
            name, fi.get("type") or "", fi.get("description") or "", fi.get("section") or ""
        )
        if http_base:
            uri = f"{http_base}/{quote(path, safe='/')}"
        elif ftp_base:
            uri = f"{ftp_base}/{quote(path, safe='/')}"
        else:
            uri = f"https://www.ebi.ac.uk/biostudies/files/{accession}/{quote(path, safe='/')}"
        artifact_id = f"ae:{accession}:{hashlib.sha1(path.encode()).hexdigest()[:16]}"
        ach = content_hash({"uri": uri, "size": fi.get("size"), "name": name})
        conn.execute(
            """
            INSERT INTO dataset_artifacts(
                artifact_id, dataset_id, sample_id, kind, role, title, filename,
                format, uri, size_bytes, access, license, source, derived_from,
                probe_status, last_checked, content_hash
            ) VALUES (?, ?, NULL, ?, ?, ?, ?, ?, ?, ?, ?, NULL, 'arrayexpress', NULL, ?, ?, ?)
            ON CONFLICT(artifact_id) DO UPDATE SET
                dataset_id=excluded.dataset_id, uri=excluded.uri, size_bytes=excluded.size_bytes,
                kind=excluded.kind, role=excluded.role, last_checked=excluded.last_checked,
                content_hash=excluded.content_hash
            """,
            (
                artifact_id, canonical_id, kind, role, fi.get("description") or name, name,
                fi.get("format") or format_from_name(name), uri, fi.get("size"), "public",
                "listed", fetched_at, ach,
            ),
        )
        cur["artifacts_created"] = cur.get("artifacts_created", 0) + 1

    if SDRF_FETCH:
        sdrf = None
        for fi in files_json:
            n = (fi.get("name") or "").lower()
            if n.endswith(".sdrf.txt") or "sdrf" in n:
                sdrf = fi
                break
        if sdrf and (sdrf.get("size") is None or sdrf["size"] <= MAX_SDRF_BYTES):
            path = sdrf.get("path") or sdrf.get("name")
            if http_base and path:
                url = f"{http_base}/{quote(path, safe='/')}"
                try:
                    body = ebi_get(url, raw=True)
                    if body and len(body) <= MAX_SDRF_BYTES:
                        samples = parse_sdrf_samples(body, canonical_id, accession)
                        for sm in samples:
                            sch = content_hash(sm["attributes"])
                            conn.execute(
                                """
                                INSERT INTO samples(
                                    sample_id, dataset_id, title, organism, tissue, disease,
                                    condition, sex, age, attributes, run_ids, content_hash, updated_at
                                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, NULL, ?, ?)
                                ON CONFLICT(sample_id) DO UPDATE SET
                                    dataset_id=excluded.dataset_id, title=excluded.title,
                                    organism=excluded.organism, tissue=excluded.tissue,
                                    disease=excluded.disease, condition=excluded.condition,
                                    sex=excluded.sex, age=excluded.age, attributes=excluded.attributes,
                                    content_hash=excluded.content_hash, updated_at=excluded.updated_at
                                """,
                                (
                                    sm["sample_id"], canonical_id, sm["title"],
                                    sm["organism"] or organism, sm["tissue"], sm["disease"],
                                    sm["condition"], sm["sex"], sm["age"],
                                    json.dumps(sm["attributes"], ensure_ascii=False), sch, fetched_at,
                                ),
                            )
                            cur["samples_created"] = cur.get("samples_created", 0) + 1
                except Exception as e:
                    log.warning("SDRF fetch/parse failed %s: %s", accession, e)

    conn.commit()


def find_main_dataset(main: sqlite3.Connection, geos: list[str], sras: list[str], prjs: list[str]) -> str | None:
    for g in geos:
        row = main.execute(
            "SELECT dataset_id FROM dataset_members WHERE accession=?", (g,)
        ).fetchone()
        if row:
            return row[0]
    for p in prjs:
        row = main.execute(
            "SELECT dataset_id FROM dataset_members WHERE accession=?", (p,)
        ).fetchone()
        if row:
            return row[0]
    for s in sras:
        row = main.execute(
            "SELECT dataset_id FROM dataset_members WHERE accession=?", (s,)
        ).fetchone()
        if row:
            return row[0]
    return None


def merge_one_accession(main: sqlite3.Connection, st: sqlite3.Connection, accession: str, cur: dict) -> None:
    src = st.execute(
        "SELECT accession, title, organism, assay, sample_count, related_geo, related_sra, files, raw_json, fetched_at FROM src_arrayexpress_entry WHERE accession=?",
        (accession,),
    ).fetchone()
    if not src:
        return
    (_acc, title, organism, assay, sample_count, related_geo, related_sra, files, raw_json, fetched_at) = src
    geos = [x for x in (related_geo or "").split(",") if x]
    sras_prjs = [x for x in (related_sra or "").split(",") if x]
    prjs = [x for x in sras_prjs if x.startswith("PRJ")]
    sras = [x for x in sras_prjs if not x.startswith("PRJ")]
    # also from raw_json
    try:
        raw = json.loads(raw_json or "{}")
        for g in raw.get("related_geo") or []:
            if g not in geos:
                geos.append(g)
        for s in raw.get("related_sra") or []:
            if s not in sras:
                sras.append(s)
        for p in raw.get("related_prj") or []:
            if p not in prjs:
                prjs.append(p)
    except Exception:
        pass

    m_geod = re.match(r"E-GEOD-(\d+)$", accession, re.I)
    if m_geod:
        gse = f"GSE{m_geod.group(1)}"
        if gse not in geos:
            geos.append(gse)

    # upsert src into main
    main.execute(
        """
        INSERT INTO src_arrayexpress_entry(
            accession, title, organism, assay, sample_count,
            related_geo, related_sra, files, raw_json, fetched_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(accession) DO UPDATE SET
            title=excluded.title, organism=excluded.organism, assay=excluded.assay,
            sample_count=excluded.sample_count, related_geo=excluded.related_geo,
            related_sra=excluded.related_sra, files=excluded.files,
            raw_json=NULL, fetched_at=excluded.fetched_at
        """,
        (_acc, title, organism, assay, sample_count, related_geo, related_sra, files, None, fetched_at),
    )

    existing_ae = main.execute(
        "SELECT dataset_id FROM dataset_members WHERE source='arrayexpress' AND accession=?",
        (accession,),
    ).fetchone()
    target = find_main_dataset(main, geos, sras, prjs)
    merged = False
    if existing_ae:
        canonical_id = existing_ae[0]
        if target and target != canonical_id:
            # prefer GEO/SRA target: re-point AE member
            canonical_id = target
            merged = True
    elif target:
        canonical_id = target
        merged = True
    else:
        canonical_id = f"ae:{accession}"
        ch = content_hash({"acc": accession, "title": title, "n": sample_count})
        prev = main.execute(
            "SELECT first_seen_at FROM datasets WHERE canonical_id=?", (canonical_id,)
        ).fetchone()
        first_seen = prev[0] if prev else (fetched_at or now_iso())
        main.execute(
            """
            INSERT INTO datasets(
                canonical_id, title, title_norm, organism, tissue, disease_tags,
                assay, platforms, n_samples, access, license, confidence,
                content_hash, first_seen_at, updated_at, probe_status, last_checked
            ) VALUES (?, ?, ?, ?, NULL, NULL, ?, NULL, ?, 'public', NULL, 'metadata_only', ?, ?, ?, 'metadata_fetched', ?)
            ON CONFLICT(canonical_id) DO UPDATE SET
                title=excluded.title, title_norm=excluded.title_norm,
                organism=COALESCE(excluded.organism, datasets.organism),
                assay=COALESCE(excluded.assay, datasets.assay),
                n_samples=COALESCE(excluded.n_samples, datasets.n_samples),
                content_hash=excluded.content_hash, updated_at=excluded.updated_at,
                last_checked=excluded.last_checked
            """,
            (
                canonical_id, title, title_norm(title), organism, assay, sample_count,
                ch, first_seen, fetched_at or now_iso(), fetched_at or now_iso(),
            ),
        )

    if merged:
        cur["merged"] = cur.get("merged", 0) + 1
        # light refresh of target dataset
        main.execute(
            """
            UPDATE datasets SET
                n_samples=COALESCE(n_samples, ?),
                assay=COALESCE(assay, ?),
                organism=COALESCE(organism, ?),
                updated_at=?
            WHERE canonical_id=?
            """,
            (sample_count, assay, organism, fetched_at or now_iso(), canonical_id),
        )

    main.execute(
        """
        INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
        VALUES (?, 'arrayexpress', ?, 'entry', ?)
        ON CONFLICT(source, accession) DO UPDATE SET
            dataset_id=excluded.dataset_id, role=excluded.role, external_url=excluded.external_url
        """,
        (
            canonical_id,
            accession,
            f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}",
        ),
    )
    for g in geos:
        main.execute(
            """
            INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
            VALUES (?, 'geo', ?, 'series', ?)
            ON CONFLICT(source, accession) DO NOTHING
            """,
            (canonical_id, g, f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}"),
        )
    for p in prjs:
        main.execute(
            """
            INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
            VALUES (?, 'sra', ?, 'bioproject', ?)
            ON CONFLICT(source, accession) DO NOTHING
            """,
            (canonical_id, p, f"https://www.ncbi.nlm.nih.gov/bioproject/{p}"),
        )
    for s in sras:
        role = "study" if re.match(r"[SED]RP\d+$", s, re.I) else "experiment"
        main.execute(
            """
            INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
            VALUES (?, 'sra', ?, ?, ?)
            ON CONFLICT(source, accession) DO NOTHING
            """,
            (canonical_id, s, role, f"https://www.ebi.ac.uk/ena/browser/view/{s}"),
        )

    # artifacts from staging provisional id
    prov = f"ae:{accession}"
    for row in st.execute(
        """SELECT artifact_id, kind, role, title, filename, format, uri, size_bytes, access,
                  license, source, derived_from, probe_status, last_checked, content_hash
           FROM dataset_artifacts WHERE dataset_id=?""",
        (prov,),
    ):
        (artifact_id, kind, role, title_a, filename, fmt, uri, size_bytes, access,
         license, source, derived_from, probe_status, last_checked, ch) = row
        main.execute(
            """
            INSERT INTO dataset_artifacts(
                artifact_id, dataset_id, sample_id, kind, role, title, filename,
                format, uri, size_bytes, access, license, source, derived_from,
                probe_status, last_checked, content_hash
            ) VALUES (?, ?, NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(artifact_id) DO UPDATE SET
                dataset_id=excluded.dataset_id, uri=excluded.uri, size_bytes=excluded.size_bytes,
                kind=excluded.kind, role=excluded.role, last_checked=excluded.last_checked,
                content_hash=excluded.content_hash
            """,
            (
                artifact_id, canonical_id, kind, role, title_a, filename, fmt, uri, size_bytes,
                access, license, source, derived_from, probe_status, last_checked, ch,
            ),
        )

    for row in st.execute(
        """SELECT sample_id, title, organism, tissue, disease, condition, sex, age,
                  attributes, run_ids, content_hash, updated_at
           FROM samples WHERE dataset_id=?""",
        (prov,),
    ):
        (sample_id, stitle, sorg, tissue, disease, condition, sex, age, attributes, run_ids, ch, updated_at) = row
        main.execute(
            """
            INSERT INTO samples(
                sample_id, dataset_id, title, organism, tissue, disease, condition, sex, age,
                attributes, run_ids, content_hash, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(sample_id) DO UPDATE SET
                dataset_id=excluded.dataset_id, title=excluded.title, organism=excluded.organism,
                tissue=excluded.tissue, disease=excluded.disease, condition=excluded.condition,
                sex=excluded.sex, age=excluded.age, attributes=excluded.attributes,
                content_hash=excluded.content_hash, updated_at=excluded.updated_at
            """,
            (
                sample_id, canonical_id, stitle, sorg, tissue, disease, condition, sex, age,
                attributes, run_ids, ch, updated_at,
            ),
        )

    cur["merged_to_main"] = cur.get("merged_to_main", 0) + 1


def merge_pending(cur: dict, limit: int = 200) -> int:
    """Merge staging accessions not yet marked into main DB.
    No-op when WRITE_MAIN_DURING_SCAN is False (scan phase).
    """
    if not WRITE_MAIN_DURING_SCAN:
        return 0
    st = connect(STAGING_DB, busy_ms=5000)
    pending = st.execute(
        """
        SELECT s.accession FROM src_arrayexpress_entry s
        LEFT JOIN ae_merge_state m ON m.accession = s.accession
        WHERE m.accession IS NULL
        ORDER BY s.fetched_at ASC
        LIMIT ?
        """,
        (limit,),
    ).fetchall()
    if not pending:
        st.close()
        return 0

    def _merge_batch(batch: list) -> int:
        main = connect(MAIN_DB, busy_ms=60000)
        n = 0
        try:
            for (acc,) in batch:
                merge_one_accession(main, st, acc, cur)
                st.execute(
                    "INSERT INTO ae_merge_state(accession, merged_at) VALUES (?, ?) "
                    "ON CONFLICT(accession) DO UPDATE SET merged_at=excluded.merged_at",
                    (acc, now_iso()),
                )
                n += 1
            main.commit()
            st.commit()
        finally:
            main.close()
        return n

    total = 0
    batch_size = 5
    consecutive_fail = 0
    for i in range(0, len(pending), batch_size):
        batch = pending[i : i + batch_size]
        try:
            total += with_lock_retry(lambda b=batch: _merge_batch(b), f"merge_batch_{batch[0][0]}")
            consecutive_fail = 0
        except Exception as e:
            consecutive_fail += 1
            log.warning("Merge batch failed (will retry later): %s", e)
            time.sleep(min(5 * consecutive_fail, 30))
            if consecutive_fail >= 3:
                break  # main still busy; leave remaining for next merge_pending call
            continue
    st.close()
    return total


def search_page(page: int, year: int | None = None) -> dict:
    params: dict[str, Any] = {
        "page": page,
        "pageSize": PAGE_SIZE,
        "sortBy": "release_date",
        "sortOrder": "ascending",
    }
    if year is not None:
        params["query"] = f"release_date:[{year}-01-01 TO {year}-12-31]"
    return ebi_get(f"{API}/arrayexpress/search", params=params)


YEAR_START = 2001
YEAR_END = 2026  # inclusive


def main() -> int:
    log.info(
        "Starting ArrayExpress SCAN-ONLY year-windows (staging=%s; EBI page offset capped ~20000)",
        STAGING_DB,
    )
    init_staging()
    st = connect(STAGING_DB, busy_ms=5000)
    cur = load_cursor(st)

    # Migrate flat-page cursor to year-window mode once.
    if cur.get("mode") != "year_windows":
        cur["mode"] = "year_windows"
        cur["year"] = YEAR_START
        cur["page"] = 1
        cur["years_done"] = []
        cur.pop("page_size", None)
        cur.pop("page_size_remap_from", None)
        cur.pop("completed_at", None)
        save_cursor(st, cur)
        log.info("Switched cursor to year_windows from year=%s page=1", YEAR_START)

    year = int(cur.get("year") or YEAR_START)
    page = int(cur.get("page") or 1)
    years_done = list(cur.get("years_done") or [])
    log.info("Resume year=%s page=%s cursor=%s", year, page, json.dumps(cur))

    since_log = 0
    LOG_EVERY = 200

    try:
        while year <= YEAR_END:
            if year in years_done:
                year += 1
                page = 1
                cur["year"] = year
                cur["page"] = page
                save_cursor(st, cur)
                continue

            first = search_page(page if page >= 1 else 1, year=year)
            if not first:
                log.error("Empty search response for year=%s page=%s", year, page)
                return 1
            total = int(first.get("totalHits") or 0)
            cur["year_total"] = total
            cur["total"] = total  # current window size
            save_cursor(st, cur)
            log.info(
                "Year %s totalHits=%s pageSize=%s resume_page=%s",
                year, total, PAGE_SIZE, page,
            )

            if total == 0:
                years_done.append(year)
                cur["years_done"] = years_done
                year += 1
                page = 1
                cur["year"] = year
                cur["page"] = page
                save_cursor(st, cur)
                continue

            current_page = page
            page_data = first

            while True:
                hits = (page_data or {}).get("hits") or []
                if not hits:
                    log.info("No more hits year=%s page=%s", year, current_page)
                    break

                for hit in hits:
                    acc = hit.get("accession")
                    if not acc or not str(acc).startswith("E-"):
                        cur["scanned"] = cur.get("scanned", 0) + 1
                        continue
                    try:
                        staging_ingest_study(st, acc, hit, cur)
                    except Exception as e:
                        cur["errors"] = cur.get("errors", 0) + 1
                        log.exception("Failed staging %s: %s", acc, e)
                        save_cursor(st, cur, error=str(e)[:500])
                    cur["scanned"] = cur.get("scanned", 0) + 1
                    since_log += 1
                    if since_log >= LOG_EVERY:
                        since_log = 0
                        save_cursor(st, cur)
                        log.info(
                            "progress year=%s page=%s/%s scanned=%s staging_new=%s skipped=%s src=%s errors=%s",
                            year,
                            current_page,
                            max(1, (total + PAGE_SIZE - 1) // PAGE_SIZE),
                            cur.get("scanned"),
                            cur.get("ingested_new"),
                            cur.get("skipped_existing"),
                            cur.get("src_upserted"),
                            cur.get("errors"),
                        )

                cur["year"] = year
                cur["page"] = current_page + 1
                save_cursor(st, cur)
                if current_page * PAGE_SIZE >= total:
                    break
                current_page += 1
                page_data = search_page(current_page, year=year)

            years_done.append(year)
            cur["years_done"] = years_done
            year += 1
            page = 1
            cur["year"] = year
            cur["page"] = page
            save_cursor(st, cur)
            log.info("Finished year window; next year=%s", year)

        cur["completed_at"] = now_iso()
        save_cursor(st, cur)
        log.info("Finished scan: %s", json.dumps(cur))
        log.info("Next: python3 arrayexpress_merge_to_main.py  # drain staging -> omics_v1.db")
    except KeyboardInterrupt:
        log.warning("Interrupted")
        cur["year"] = year
        cur["page"] = current_page if "current_page" in locals() else page
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
