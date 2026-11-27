#!/usr/bin/env python3
"""Parquet-native, serverless incremental sync & DuckDB dedup merger.

Eliminates the need to store/maintain the 15+ GB SQLite main DB on the server:
1. Reads tiered schedule from `pipeline/sync_policy.json` and per-source watermarks
   from `serving/.sync_state.json`.
2. Fetches only records newer than `last_synced_at` (with configurable overlap_days)
   into a tiny temporary delta SQLite (deleted after merge), OR ingests an external
   `--delta-sqlite` file.
3. Performs two-stage deduplication directly inside DuckDB against `serving/current/`:
   - L1 cross-repository dedup via `dataset_members.parquet` (source, accession / accession)
     and DOI-based `publications.parquet` dedup.
   - Primary-key upsert via `QUALIFY ROW_NUMBER()` on compact tables and ONLY the
     touched `source=*` partitions of `samples` and `dataset_artifacts`.
   - Unchanged partitions are zero-copy hardlinked from the previous release.
4. Atomically publishes `serving/releases/v_YYYYMMDD_HHMMSS` -> `serving/current`
   and advances `serving/.sync_state.json` only after validation succeeds.

Usage:
  python3 pipeline/incremental_parquet_sync.py --status
  python3 pipeline/incremental_parquet_sync.py --dry-run
  python3 pipeline/incremental_parquet_sync.py --source geo --force --max-items 20
  python3 pipeline/incremental_parquet_sync.py --delta-sqlite /tmp/delta.db --source custom
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
import shutil
import sqlite3
import sys
import tempfile
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

try:
    from .paths import CODE_ROOT, POLICY_PATH, SCHEMA_PATH, serving_dir
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from paths import CODE_ROOT, POLICY_PATH, SCHEMA_PATH, serving_dir  # type: ignore

ROOT = CODE_ROOT
SERVING_DIR = serving_dir()
RELEASES_DIR = SERVING_DIR / "releases"
CURRENT_LINK = SERVING_DIR / "current"
STATE_FILE = SERVING_DIR / ".sync_state.json"
TZ = timezone(timedelta(hours=8))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s [inc-sync] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("incremental_parquet_sync")


def now_dt() -> datetime:
    return datetime.now(TZ)


def now_iso() -> str:
    return now_dt().strftime("%Y-%m-%dT%H:%M:%S%z")


def now_tag() -> str:
    return now_dt().strftime("v_%Y%m%d_%H%M%S")


def chash(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]


def tnorm(text: str | None) -> str:
    return re.sub(r"\s+", " ", (text or "").lower().strip())


def load_policy() -> dict[str, Any]:
    if not POLICY_PATH.exists():
        raise FileNotFoundError(f"Policy file not found: {POLICY_PATH}")
    return json.loads(POLICY_PATH.read_text())


def load_state() -> dict[str, Any]:
    if STATE_FILE.exists():
        return json.loads(STATE_FILE.read_text())
    # Bootstrap watermarks from serving/current/_manifest.json if available
    manifest = CURRENT_LINK / "_manifest.json"
    boot_ts = now_iso()
    if manifest.exists():
        try:
            boot_ts = json.loads(manifest.read_text()).get("created_at") or boot_ts
        except Exception:
            pass
    return {"bootstrapped_at": boot_ts, "sources": {}}


def save_state(state: dict[str, Any]) -> None:
    SERVING_DIR.mkdir(parents=True, exist_ok=True)
    tmp = STATE_FILE.with_suffix(".tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True))
    os.replace(tmp, STATE_FILE)


def parse_iso(ts: str | None) -> datetime | None:
    if not ts:
        return None
    try:
        return datetime.fromisoformat(ts.replace("Z", "+00:00"))
    except ValueError:
        try:
            return datetime.strptime(ts, "%Y-%m-%dT%H:%M:%S%z")
        except ValueError:
            return None


def is_due(cfg: dict[str, Any], src_state: dict[str, Any], force: bool) -> tuple[bool, str]:
    if force:
        return True, "forced"
    if not cfg.get("enabled", False):
        return False, "disabled"
    interval_h = float(cfg.get("interval_hours") or 0)
    if interval_h <= 0:
        return False, "manual_only"
    last_ts = parse_iso(src_state.get("last_synced_at"))
    if last_ts is None:
        return True, "never_synced"
    elapsed_h = (now_dt() - last_ts.astimezone(TZ)).total_seconds() / 3600.0
    if elapsed_h >= interval_h:
        return True, f"elapsed={elapsed_h:.1f}h>={interval_h}h"
    return False, f"elapsed={elapsed_h:.1f}h<{interval_h}h"


def init_delta_db(path: Path) -> sqlite3.Connection:
    if path.exists():
        path.unlink()
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys = OFF")
    conn.executescript(SCHEMA_PATH.read_text())
    conn.commit()
    return conn


def _http_get_json(url: str, params: dict | None = None, timeout: int = 30) -> Any:
    full_url = f"{url}?{urlencode(params)}" if params else url
    for attempt in range(1, 5):
        try:
            req = Request(full_url, headers={"User-Agent": "omics-sniffer/1.0 (incremental-sync)"})
            with urlopen(req, timeout=timeout) as resp:
                if resp.status in (429,) or resp.status >= 500:
                    time.sleep(min(2 ** attempt, 16))
                    continue
                return json.loads(resp.read().decode("utf-8"))
        except Exception as e:
            if attempt == 4:
                raise
            log.warning("HTTP retry %d/4 %s: %s", attempt, url, e)
            time.sleep(min(2 ** attempt, 16))
    return None


# ---------------------------------------------------------------------------
# Lightweight incremental fetchers for Tier-1 daily sources
# ---------------------------------------------------------------------------

def fetch_geo_delta(conn: sqlite3.Connection, since_dt: datetime, max_items: int | None) -> int:
    """Fetch GEO Series (GSE) published/updated since `since_dt`."""
    start_str = since_dt.strftime("%Y/%m/%d")
    end_str = (now_dt() + timedelta(days=1)).strftime("%Y/%m/%d")
    term = f'gse[Entry Type] AND ("{start_str}"[PDAT] : "{end_str}"[PDAT])'
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    common = {"db": "gds", "tool": "omics-sniffer", "email": "omics-sniffer@localhost", "retmode": "json"}

    meta = _http_get_json(f"{base}/esearch.fcgi", {**common, "term": term, "retmax": 0})
    total = int((meta or {}).get("esearchresult", {}).get("count") or 0)
    cap = min(total, max_items) if max_items else total
    log.info("GEO date window %s..%s: total=%d (cap=%d)", start_str, end_str, total, cap)
    if cap == 0:
        return 0

    fetched = 0
    retstart = 0
    batch = 100
    ts = now_iso()

    while retstart < cap:
        take = min(batch, cap - retstart)
        time.sleep(0.35)
        sr = _http_get_json(f"{base}/esearch.fcgi", {**common, "term": term, "retstart": retstart, "retmax": take})
        ids = (sr or {}).get("esearchresult", {}).get("idlist") or []
        if not ids:
            break
        time.sleep(0.35)
        sm = _http_get_json(f"{base}/esummary.fcgi", {**common, "id": ",".join(ids)})
        res = (sm or {}).get("result") or {}
        for uid in res.get("uids") or ids:
            doc = res.get(str(uid))
            if not isinstance(doc, dict) or doc.get("entrytype") != "GSE":
                continue
            acc = doc.get("accession") or (f"GSE{doc.get('gse')}" if str(doc.get("gse") or "").isdigit() else "")
            if not acc.startswith("GSE"):
                continue
            cid = f"geo:{acc}"
            title = doc.get("title") or ""
            organism = doc.get("taxon") or doc.get("samplestaxa") or None
            gpl = doc.get("gpl") or ""
            platforms = [p.strip() if p.strip().upper().startswith("GPL") else f"GPL{p.strip()}" for p in str(gpl).split(";") if p.strip()]
            n_samples = int(doc.get("n_samples") or 0)
            assay = doc.get("gdstype") or None
            ftplink = (doc.get("ftplink") or "").rstrip("/") + ("/" if doc.get("ftplink") else "")
            ch = chash({"acc": acc, "title": title, "n": n_samples, "assay": assay})

            conn.execute(
                """
                INSERT OR REPLACE INTO datasets(
                    canonical_id, title, title_norm, organism, tissue, disease_tags,
                    assay, platforms, n_samples, access, license, confidence,
                    content_hash, first_seen_at, updated_at, probe_status, last_checked
                ) VALUES (?, ?, ?, ?, NULL, NULL, ?, ?, ?, 'public', NULL, 'metadata_only', ?, ?, ?, 'metadata_fetched', ?)
                """,
                (cid, title, tnorm(title), organism, assay, json.dumps(platforms), n_samples, ch, ts, ts, ts),
            )
            conn.execute(
                "INSERT OR IGNORE INTO dataset_members(dataset_id, source, accession, role, external_url) VALUES (?, 'geo', ?, 'series', ?)",
                (cid, acc, f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"),
            )
            if doc.get("bioproject"):
                bp = str(doc["bioproject"]).strip()
                conn.execute(
                    "INSERT OR IGNORE INTO dataset_members(dataset_id, source, accession, role, external_url) VALUES (?, 'ncbi', ?, 'bioproject', ?)",
                    (cid, bp, f"https://www.ncbi.nlm.nih.gov/bioproject/{bp}"),
                )
            for s in doc.get("samples") or []:
                if isinstance(s, dict) and str(s.get("accession") or "").startswith("GSM"):
                    gsm = s["accession"]
                    stitle = s.get("title") or ""
                    conn.execute(
                        """
                        INSERT OR REPLACE INTO samples(
                            sample_id, dataset_id, title, organism, tissue, disease,
                            condition, sex, age, attributes, run_ids, content_hash, updated_at
                        ) VALUES (?, ?, ?, ?, NULL, NULL, NULL, NULL, NULL, ?, NULL, ?, ?)
                        """,
                        (f"geo:{gsm}", cid, stitle, organism, json.dumps({"series": acc}), chash({"gsm": gsm, "title": stitle}), ts),
                    )
            if ftplink:
                conn.execute(
                    """
                    INSERT OR REPLACE INTO dataset_artifacts(
                        artifact_id, dataset_id, sample_id, kind, role, title, filename,
                        format, uri, size_bytes, access, license, source, derived_from,
                        probe_status, last_checked, content_hash
                    ) VALUES (?, ?, NULL, 'raw', 'ftp_root', ?, NULL, NULL, ?, NULL, 'public', NULL, 'geo', NULL, 'uri_only', ?, ?)
                    """,
                    (f"geo:{acc}:ftp_root", cid, f"{acc} GEO FTP root", ftplink, ts, chash({"uri": ftplink})),
                )
            for pmid in doc.get("pubmedids") or []:
                pmid_s = str(pmid).strip()
                if pmid_s.isdigit():
                    conn.execute(
                        "INSERT OR IGNORE INTO publications(pmid, url, fetched_at) VALUES (?, ?, ?)",
                        (pmid_s, f"https://pubmed.ncbi.nlm.nih.gov/{pmid_s}/", ts),
                    )
                    conn.execute(
                        "INSERT OR IGNORE INTO dataset_publications(dataset_id, pmid, relation, source) VALUES (?, ?, 'primary', 'geo_link')",
                        (cid, pmid_s),
                    )
            fetched += 1
        conn.commit()
        retstart += len(ids)
    return fetched


def fetch_arrayexpress_delta(conn: sqlite3.Connection, since_dt: datetime, max_items: int | None) -> int:
    """Fetch ArrayExpress studies released since `since_dt` via BioStudies search."""
    start_str = since_dt.strftime("%Y-%m-%d")
    end_str = (now_dt() + timedelta(days=1)).strftime("%Y-%m-%d")
    api = "https://www.ebi.ac.uk/biostudies/api/v1/arrayexpress/search"
    query = f"release_date:[{start_str} TO {end_str}]"

    first = _http_get_json(api, {"query": query, "page": 1, "pageSize": 20, "sortBy": "release_date", "sortOrder": "descending"})
    total = int((first or {}).get("totalHits") or 0)
    cap = min(total, max_items) if max_items else total
    log.info("ArrayExpress window %s..%s: total=%d (cap=%d)", start_str, end_str, total, cap)
    if cap == 0:
        return 0

    fetched = 0
    page = 1
    page_data = first
    ts = now_iso()

    while fetched < cap and page_data:
        hits = page_data.get("hits") or []
        if not hits:
            break
        for h in hits:
            if fetched >= cap:
                break
            acc = str(h.get("accession") or "").strip()
            if not acc.startswith("E-"):
                continue
            cid = f"ae:{acc}"
            title = h.get("title") or acc
            ch = chash({"acc": acc, "title": title})
            conn.execute(
                """
                INSERT OR REPLACE INTO datasets(
                    canonical_id, title, title_norm, organism, tissue, disease_tags,
                    assay, platforms, n_samples, access, license, confidence,
                    content_hash, first_seen_at, updated_at, probe_status, last_checked
                ) VALUES (?, ?, ?, NULL, NULL, NULL, NULL, NULL, NULL, ?, NULL, 'metadata_only', ?, ?, ?, 'metadata_fetched', ?)
                """,
                (cid, title, tnorm(title), "public" if h.get("isPublic", True) else "unknown", ch, ts, ts, ts),
            )
            conn.execute(
                "INSERT OR IGNORE INTO dataset_members(dataset_id, source, accession, role, external_url) VALUES (?, 'arrayexpress', ?, 'entry', ?)",
                (cid, acc, f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{acc}"),
            )
            m_geod = re.match(r"E-GEOD-(\d+)$", acc, re.I)
            if m_geod:
                gse = f"GSE{m_geod.group(1)}"
                conn.execute(
                    "INSERT OR IGNORE INTO dataset_members(dataset_id, source, accession, role, external_url) VALUES (?, 'geo', ?, 'series_link', ?)",
                    (cid, gse, f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"),
                )
            fetched += 1
        conn.commit()
        page += 1
        if fetched >= cap:
            break
        time.sleep(0.25)
        page_data = _http_get_json(api, {"query": query, "page": page, "pageSize": 20, "sortBy": "release_date", "sortOrder": "descending"})
    return fetched


def fetch_pride_delta(conn: sqlite3.Connection, since_dt: datetime, max_items: int | None) -> int:
    """Fetch PRIDE Archive projects sorted descending by publicationDate with early stop."""
    since_str = since_dt.strftime("%Y-%m-%d")
    api = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects"
    fetched = 0
    page = 0
    ts = now_iso()

    while True:
        data = _http_get_json(api, {"page": page, "pageSize": 50, "sortConditions": "publicationDate", "sortDirection": "DESC"})
        if not isinstance(data, list) or not data:
            break
        stop = False
        for p in data:
            if max_items and fetched >= max_items:
                stop = True
                break
            pub_date = str(p.get("publicationDate") or "")[:10]
            if pub_date and pub_date < since_str:
                stop = True
                break
            acc = str(p.get("accession") or "").strip().upper()
            if not acc.startswith("PXD"):
                continue
            cid = f"pride:{acc}"
            title = p.get("title") or acc
            organisms = "; ".join(str(x.get("name") or x) for x in (p.get("organisms") or []) if x) or None
            diseases = [str(x.get("name") or x) for x in (p.get("diseases") or []) if x]
            ch = chash({"acc": acc, "title": title, "pub": pub_date})
            conn.execute(
                """
                INSERT OR REPLACE INTO datasets(
                    canonical_id, title, title_norm, organism, tissue, disease_tags,
                    assay, platforms, n_samples, access, license, confidence,
                    content_hash, first_seen_at, updated_at, probe_status, last_checked
                ) VALUES (?, ?, ?, ?, NULL, ?, 'proteomics', NULL, NULL, 'public', ?, 'metadata_only', ?, ?, ?, 'metadata_fetched', ?)
                """,
                (cid, title, tnorm(title), organisms, json.dumps(diseases) if diseases else None, p.get("license"), ch, ts, ts, ts),
            )
            conn.execute(
                "INSERT OR IGNORE INTO dataset_members(dataset_id, source, accession, role, external_url) VALUES (?, 'pride', ?, 'project', ?)",
                (cid, acc, f"https://www.ebi.ac.uk/pride/archive/projects/{acc}"),
            )
            fetched += 1
        conn.commit()
        if stop:
            break
        page += 1
        time.sleep(0.3)
    log.info("PRIDE reverse early-stop since %s: fetched=%d", since_str, fetched)
    return fetched


def fetch_cellxgene_delta(conn: sqlite3.Connection, since_dt: datetime, max_items: int | None) -> int:
    """Fetch CZ CELLxGENE Curation datasets published or revised since `since_dt`."""
    since_str = since_dt.strftime("%Y-%m-%d")
    api = "https://api.cellxgene.cziscience.com/curation/v1/datasets"
    data = _http_get_json(api, timeout=60)
    if not isinstance(data, list):
        return 0
    fetched = 0
    ts = now_iso()
    for d in data:
        if max_items and fetched >= max_items:
            break
        if not isinstance(d, dict):
            continue
        pub_at = str(d.get("published_at") or "")[:10]
        rev_at = str(d.get("revised_at") or "")[:10]
        latest_dt = max(pub_at, rev_at)
        if latest_dt and latest_dt < since_str:
            continue
        ds_uuid = str(d.get("dataset_id") or "").strip()
        if not ds_uuid:
            continue
        cid = f"cellxgene:{ds_uuid}"
        title = d.get("title") or ds_uuid
        orgs = "; ".join(str(x.get("label") or "") for x in (d.get("organism") or []) if isinstance(x, dict) and x.get("label")) or None
        tissues = "; ".join(str(x.get("label") or "") for x in (d.get("tissue") or []) if isinstance(x, dict) and x.get("label")) or None
        diseases = [str(x.get("label")) for x in (d.get("disease") or []) if isinstance(x, dict) and x.get("label")]
        assays = "; ".join(str(x.get("label") or "") for x in (d.get("assay") or []) if isinstance(x, dict) and x.get("label")) or "scRNA-seq"
        n_cells = int(d.get("cell_count") or 0) or None
        ch = chash({"id": ds_uuid, "title": title, "rev": rev_at or pub_at, "n": n_cells})
        conn.execute(
            """
            INSERT OR REPLACE INTO datasets(
                canonical_id, title, title_norm, organism, tissue, disease_tags,
                assay, platforms, n_samples, access, license, confidence,
                content_hash, first_seen_at, updated_at, probe_status, last_checked
            ) VALUES (?, ?, ?, ?, ?, ?, ?, NULL, ?, 'public', 'CC-BY-4.0', 'metadata_only', ?, ?, ?, 'metadata_fetched', ?)
            """,
            (cid, title, tnorm(title), orgs, tissues, json.dumps(diseases) if diseases else None, assays, n_cells, ch, ts, ts, ts),
        )
        conn.execute(
            "INSERT OR IGNORE INTO dataset_members(dataset_id, source, accession, role, external_url) VALUES (?, 'cellxgene', ?, 'dataset', ?)",
            (cid, ds_uuid, f"https://cellxgene.cziscience.com/e/{ds_uuid}.cxg/"),
        )
        fetched += 1
    conn.commit()
    log.info("CELLxGENE window since %s: fetched=%d", since_str, fetched)
    return fetched


def fetch_metabolights_delta(conn: sqlite3.Connection, since_dt: datetime, max_items: int | None) -> int:
    """Fetch MetaboLights studies released since `since_dt` via EBI Search."""
    start_str = since_dt.strftime("%Y%m%d")
    end_str = (now_dt() + timedelta(days=1)).strftime("%Y%m%d")
    api = "https://www.ebi.ac.uk/ebisearch/ws/rest/metabolights"
    query = f"publication_date:[{start_str} TO {end_str}]"
    fetched = 0
    start_offset = 0
    size = 100
    ts = now_iso()

    while True:
        res = _http_get_json(
            api,
            {
                "query": query,
                "format": "json",
                "size": size,
                "start": start_offset,
                "fields": "name,description,omics_type,Organism,Organism_Part",
            },
        )
        entries = (res or {}).get("entries") or []
        hit_count = int((res or {}).get("hitCount") or 0)
        cap = min(hit_count, max_items) if max_items else hit_count
        if not entries or fetched >= cap:
            break
        for e in entries:
            if fetched >= cap:
                break
            acc = str(e.get("id") or "").strip().upper()
            if not acc.startswith("MTBLS"):
                continue
            fields = e.get("fields") or {}
            title = ((fields.get("name") or [acc])[0]) if isinstance(fields.get("name"), list) else acc
            org = "; ".join(fields.get("Organism") or []) or None
            tissue = "; ".join(fields.get("Organism_Part") or []) or None
            cid = f"mtbls:{acc}"
            ch = chash({"acc": acc, "title": title, "org": org})
            conn.execute(
                """
                INSERT OR REPLACE INTO datasets(
                    canonical_id, title, title_norm, organism, tissue, disease_tags,
                    assay, platforms, n_samples, access, license, confidence,
                    content_hash, first_seen_at, updated_at, probe_status, last_checked
                ) VALUES (?, ?, ?, ?, ?, NULL, 'metabolomics', NULL, NULL, 'public', NULL, 'metadata_only', ?, ?, ?, 'metadata_fetched', ?)
                """,
                (cid, title, tnorm(title), org, tissue, ch, ts, ts, ts),
            )
            conn.execute(
                "INSERT OR IGNORE INTO dataset_members(dataset_id, source, accession, role, external_url) VALUES (?, 'metabolights', ?, 'study', ?)",
                (cid, acc, f"https://www.ebi.ac.uk/metabolights/{acc}"),
            )
            fetched += 1
        conn.commit()
        start_offset += len(entries)
        if start_offset >= cap:
            break
        time.sleep(0.25)
    log.info("MetaboLights window %s..%s: fetched=%d", start_str, end_str, fetched)
    return fetched


FETCHERS = {
    "geo": fetch_geo_delta,
    "arrayexpress": fetch_arrayexpress_delta,
    "pride": fetch_pride_delta,
    "cellxgene": fetch_cellxgene_delta,
    "metabolights": fetch_metabolights_delta,
}


# ---------------------------------------------------------------------------
# DuckDB two-stage dedup merge (Parquet -> Parquet, atomic swap)
# ---------------------------------------------------------------------------

def link_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def atomic_symlink_swap(target_dir: Path, link_path: Path) -> None:
    rel_target = os.path.relpath(target_dir, link_path.parent)
    tmp_link = link_path.parent / f".{link_path.name}.tmp.{os.getpid()}"
    if tmp_link.exists() or tmp_link.is_symlink():
        tmp_link.unlink()
    os.symlink(rel_target, tmp_link)
    os.replace(tmp_link, link_path)


def cleanup_old_releases(keep: int) -> None:
    if not RELEASES_DIR.exists():
        return
    current_target = CURRENT_LINK.resolve() if CURRENT_LINK.exists() else None
    releases = sorted(
        [d for d in RELEASES_DIR.iterdir() if d.is_dir() and d.name.startswith("v_")],
        key=lambda d: d.name,
        reverse=True,
    )
    for old in releases[keep:]:
        if current_target and old.resolve() == current_target:
            continue
        log.info("Removing old release: %s", old.name)
        shutil.rmtree(old, ignore_errors=True)


def merge_delta_into_parquet(delta_db: Path, memory_limit: str, keep_releases: int) -> dict[str, Any]:
    import duckdb  # type: ignore

    if not CURRENT_LINK.exists():
        raise FileNotFoundError(f"Baseline Parquet release not found at {CURRENT_LINK}. Run export_parquet.py once first.")

    prev_release = CURRENT_LINK.resolve()
    release_dir = RELEASES_DIR / now_tag()
    release_dir.mkdir(parents=True, exist_ok=True)
    spill_dir = Path(tempfile.mkdtemp(prefix="duckdb_spill_"))

    con = duckdb.connect(database=":memory:")
    try:
        con.execute(f"SET memory_limit = '{memory_limit}';")
        con.execute(f"SET temp_directory = '{spill_dir}';")
        con.execute("INSTALL sqlite; LOAD sqlite;")
        con.execute(f"ATTACH '{delta_db}' AS delta (TYPE SQLITE, READ_ONLY);")

        delta_ds_cnt = con.execute("SELECT COUNT(*) FROM delta.datasets").fetchone()[0]
        if delta_ds_cnt == 0:
            shutil.rmtree(release_dir, ignore_errors=True)
            return {"status": "no_changes", "delta_datasets": 0}

        # ------------------------------------------------------------------
        # Stage A: L1 cross-repo dedup via dataset_members.parquet & DOI map
        # ------------------------------------------------------------------
        mem_pq = str(prev_release / "dataset_members.parquet")
        pub_pq = str(prev_release / "publications.parquet")

        con.execute(
            f"""
            CREATE TEMP TABLE id_map AS
            WITH exact_hit AS (
                SELECT dm.dataset_id AS st_id, cur.dataset_id AS target_id, 1 AS prio
                FROM delta.dataset_members dm
                JOIN read_parquet('{mem_pq}') cur
                  ON cur.source = dm.source AND cur.accession = dm.accession
            ),
            acc_hit AS (
                SELECT dm.dataset_id AS st_id, cur.dataset_id AS target_id, 2 AS prio
                FROM delta.dataset_members dm
                JOIN read_parquet('{mem_pq}') cur
                  ON cur.accession = dm.accession
            ),
            all_hits AS (
                SELECT * FROM exact_hit
                UNION ALL
                SELECT * FROM acc_hit
            ),
            best_hit AS (
                SELECT st_id, target_id
                FROM all_hits
                QUALIFY ROW_NUMBER() OVER (PARTITION BY st_id ORDER BY prio ASC, target_id ASC) = 1
            )
            SELECT d.canonical_id AS st_id, COALESCE(b.target_id, d.canonical_id) AS target_id
            FROM delta.datasets d
            LEFT JOIN best_hit b ON b.st_id = d.canonical_id;
            """
        )

        l1_deduped = con.execute("SELECT COUNT(*) FROM id_map WHERE st_id != target_id").fetchone()[0]

        con.execute(
            f"""
            CREATE TEMP TABLE pmid_map AS
            SELECT dp.pmid AS st_pmid, COALESCE(cur.pmid, dp.pmid) AS target_pmid
            FROM delta.publications dp
            LEFT JOIN read_parquet('{pub_pq}') cur
              ON dp.doi IS NOT NULL AND cur.doi = dp.doi;
            """
        )

        # ------------------------------------------------------------------
        # Stage B: Compact tables (datasets, publications, dataset_publications, dataset_members)
        # ------------------------------------------------------------------
        ds_pq = str(prev_release / "datasets.parquet")
        prev_ds_cnt = con.execute(f"SELECT COUNT(*) FROM read_parquet('{ds_pq}')").fetchone()[0]

        con.execute(
            f"""
            COPY (
                WITH delta_norm AS (
                    SELECT m.target_id AS canonical_id, d.title, d.title_norm, d.organism, d.tissue,
                           d.disease_tags, d.assay, d.platforms, d.n_samples, d.access, d.license,
                           d.confidence, d.content_hash, d.first_seen_at, d.updated_at, d.probe_status, d.last_checked,
                           CASE WHEN m.st_id = m.target_id THEN 1 ELSE 3 END AS _prio
                    FROM delta.datasets d
                    JOIN id_map m ON m.st_id = d.canonical_id
                ),
                merged AS (
                    SELECT *, 2 AS _prio FROM read_parquet('{ds_pq}')
                    UNION ALL
                    SELECT * FROM delta_norm
                )
                SELECT canonical_id, title, title_norm, organism, tissue, disease_tags,
                       assay, platforms, n_samples, access, license, confidence,
                       content_hash, first_seen_at, updated_at, probe_status, last_checked
                FROM merged
                QUALIFY ROW_NUMBER() OVER (PARTITION BY canonical_id ORDER BY _prio ASC, COALESCE(updated_at, '') DESC) = 1
            ) TO '{release_dir / "datasets.parquet"}' (FORMAT PARQUET, COMPRESSION 'ZSTD', ROW_GROUP_SIZE 100000);
            """
        )

        delta_pub_cnt = con.execute("SELECT COUNT(*) FROM delta.publications").fetchone()[0]
        if delta_pub_cnt > 0:
            con.execute(
                f"""
                COPY (
                    WITH delta_norm AS (
                        SELECT pm.target_pmid AS pmid, p.doi, p.title, p.authors, p.journal,
                               p.year, p.abstract, p.url, p.fetched_at, 1 AS _prio
                        FROM delta.publications p
                        JOIN pmid_map pm ON pm.st_pmid = p.pmid
                    ),
                    merged AS (
                        SELECT *, 2 AS _prio FROM read_parquet('{pub_pq}')
                        UNION ALL
                        SELECT * FROM delta_norm
                    )
                    SELECT pmid, doi, title, authors, journal, year, abstract, url, fetched_at
                    FROM merged
                    QUALIFY ROW_NUMBER() OVER (PARTITION BY pmid ORDER BY _prio ASC, COALESCE(fetched_at, '') DESC) = 1
                ) TO '{release_dir / "publications.parquet"}' (FORMAT PARQUET, COMPRESSION 'ZSTD', ROW_GROUP_SIZE 100000);
                """
            )
        else:
            link_or_copy(prev_release / "publications.parquet", release_dir / "publications.parquet")

        dp_pq = str(prev_release / "dataset_publications.parquet")
        delta_dp_cnt = con.execute("SELECT COUNT(*) FROM delta.dataset_publications").fetchone()[0]
        if delta_dp_cnt > 0:
            con.execute(
                f"""
                COPY (
                    WITH merged AS (
                        SELECT dataset_id, pmid, relation, source, 2 AS _prio FROM read_parquet('{dp_pq}')
                        UNION ALL
                        SELECT m.target_id AS dataset_id, COALESCE(pm.target_pmid, dp.pmid) AS pmid,
                               dp.relation, dp.source, 1 AS _prio
                        FROM delta.dataset_publications dp
                        JOIN id_map m ON m.st_id = dp.dataset_id
                        LEFT JOIN pmid_map pm ON pm.st_pmid = dp.pmid
                    )
                    SELECT dataset_id, pmid, relation, source
                    FROM merged
                    QUALIFY ROW_NUMBER() OVER (PARTITION BY dataset_id, pmid ORDER BY _prio ASC) = 1
                ) TO '{release_dir / "dataset_publications.parquet"}' (FORMAT PARQUET, COMPRESSION 'ZSTD', ROW_GROUP_SIZE 100000);
                """
            )
        else:
            link_or_copy(prev_release / "dataset_publications.parquet", release_dir / "dataset_publications.parquet")

        delta_mem_cnt = con.execute("SELECT COUNT(*) FROM delta.dataset_members").fetchone()[0]
        if delta_mem_cnt > 0:
            con.execute(
                f"""
                COPY (
                    WITH max_id AS (SELECT COALESCE(MAX(id), 0) AS m FROM read_parquet('{mem_pq}')),
                    delta_norm AS (
                        SELECT (SELECT m FROM max_id) + ROW_NUMBER() OVER () AS id,
                               m.target_id AS dataset_id, dm.source, dm.accession, dm.role, dm.external_url, 2 AS _prio
                        FROM delta.dataset_members dm
                        JOIN id_map m ON m.st_id = dm.dataset_id
                    ),
                    merged AS (
                        SELECT id, dataset_id, source, accession, role, external_url, 1 AS _prio FROM read_parquet('{mem_pq}')
                        UNION ALL
                        SELECT * FROM delta_norm
                    )
                    SELECT id, dataset_id, source, accession, role, external_url
                    FROM merged
                    QUALIFY ROW_NUMBER() OVER (PARTITION BY source, accession ORDER BY _prio ASC) = 1
                ) TO '{release_dir / "dataset_members.parquet"}' (FORMAT PARQUET, COMPRESSION 'ZSTD', ROW_GROUP_SIZE 100000);
                """
            )
        else:
            link_or_copy(prev_release / "dataset_members.parquet", release_dir / "dataset_members.parquet")

        # ------------------------------------------------------------------
        # Stage C: Partitioned tables (samples & dataset_artifacts)
        # Only rewrite touched source=* partitions; hardlink all others
        # ------------------------------------------------------------------
        # 1. samples
        touched_sample_sources = {
            r[0]
            for r in con.execute(
                """
                SELECT DISTINCT
                    CASE WHEN instr(sample_id, ':') > 0
                         THEN split_part(sample_id, ':', 1)
                         ELSE 'geo'
                    END AS src
                FROM delta.samples
                """
            ).fetchall()
        }
        prev_samples_dir = prev_release / "samples"
        new_samples_dir = release_dir / "samples"
        new_samples_dir.mkdir(parents=True, exist_ok=True)
        existing_sample_sources = {
            d.name.split("=", 1)[1] for d in prev_samples_dir.glob("source=*") if d.is_dir()
        }
        for src in sorted(existing_sample_sources | touched_sample_sources):
            dst = new_samples_dir / f"source={src}" / "data.parquet"
            old_file = prev_samples_dir / f"source={src}" / "data.parquet"
            if src not in touched_sample_sources:
                link_or_copy(old_file, dst)
                continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            cond = "instr(s.sample_id, ':') = 0" if src == "geo" else f"split_part(s.sample_id, ':', 1) = '{src}'"
            old_cte = (
                f"SELECT sample_id, dataset_id, title, organism, tissue, disease, condition, sex, age, attributes, run_ids, content_hash, updated_at, 2 AS _prio FROM read_parquet('{old_file}') UNION ALL"
                if old_file.exists()
                else ""
            )
            log.info("Rewriting samples/source=%s ...", src)
            con.execute(
                f"""
                COPY (
                    WITH merged AS (
                        {old_cte}
                        SELECT s.sample_id, m.target_id AS dataset_id, s.title, s.organism, s.tissue,
                               s.disease, s.condition, s.sex, s.age, s.attributes, s.run_ids,
                               s.content_hash, s.updated_at, 1 AS _prio
                        FROM delta.samples s
                        JOIN id_map m ON m.st_id = s.dataset_id
                        WHERE {cond}
                    )
                    SELECT sample_id, dataset_id, title, organism, tissue, disease, condition,
                           sex, age, attributes, run_ids, content_hash, updated_at
                    FROM merged
                    QUALIFY ROW_NUMBER() OVER (PARTITION BY sample_id ORDER BY _prio ASC, COALESCE(updated_at, '') DESC) = 1
                ) TO '{dst}' (FORMAT PARQUET, COMPRESSION 'ZSTD', ROW_GROUP_SIZE 100000);
                """
            )

        # 2. dataset_artifacts
        touched_art_sources = {
            r[0]
            for r in con.execute(
                "SELECT DISTINCT COALESCE(NULLIF(source, ''), 'unknown') FROM delta.dataset_artifacts"
            ).fetchall()
        }
        prev_art_dir = prev_release / "dataset_artifacts"
        new_art_dir = release_dir / "dataset_artifacts"
        new_art_dir.mkdir(parents=True, exist_ok=True)
        existing_art_sources = {
            d.name.split("=", 1)[1] for d in prev_art_dir.glob("source=*") if d.is_dir()
        }
        for src in sorted(existing_art_sources | touched_art_sources):
            dst = new_art_dir / f"source={src}" / "data.parquet"
            old_file = prev_art_dir / f"source={src}" / "data.parquet"
            if src not in touched_art_sources:
                link_or_copy(old_file, dst)
                continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            cond = "(a.source IS NULL OR a.source = '')" if src == "unknown" else f"a.source = '{src}'"
            old_cte = (
                f"SELECT artifact_id, dataset_id, sample_id, kind, role, title, filename, format, uri, size_bytes, access, license, source, derived_from, probe_status, last_checked, content_hash, 2 AS _prio FROM read_parquet('{old_file}') UNION ALL"
                if old_file.exists()
                else ""
            )
            log.info("Rewriting dataset_artifacts/source=%s ...", src)
            con.execute(
                f"""
                COPY (
                    WITH merged AS (
                        {old_cte}
                        SELECT a.artifact_id, m.target_id AS dataset_id, a.sample_id, a.kind, a.role,
                               a.title, a.filename, a.format, a.uri, a.size_bytes, a.access, a.license,
                               a.source, a.derived_from, a.probe_status, a.last_checked, a.content_hash, 1 AS _prio
                        FROM delta.dataset_artifacts a
                        JOIN id_map m ON m.st_id = a.dataset_id
                        WHERE {cond}
                    )
                    SELECT artifact_id, dataset_id, sample_id, kind, role, title, filename,
                           format, uri, size_bytes, access, license, source, derived_from,
                           probe_status, last_checked, content_hash
                    FROM merged
                    QUALIFY ROW_NUMBER() OVER (PARTITION BY artifact_id ORDER BY _prio ASC, COALESCE(last_checked, '') DESC) = 1
                ) TO '{dst}' (FORMAT PARQUET, COMPRESSION 'ZSTD', ROW_GROUP_SIZE 100000);
                """
            )

        # ------------------------------------------------------------------
        # Stage D: Validate + Atomic Swap
        # ------------------------------------------------------------------
        new_ds_pq = release_dir / "datasets.parquet"
        new_pub_pq = release_dir / "publications.parquet"
        new_ds_cnt = con.execute(f"SELECT COUNT(*) FROM read_parquet('{new_ds_pq}')").fetchone()[0]
        new_pub_cnt = con.execute(f"SELECT COUNT(*) FROM read_parquet('{new_pub_pq}')").fetchone()[0]
        if new_ds_cnt < prev_ds_cnt:
            raise RuntimeError(f"Validation failed: datasets count regressed {new_ds_cnt} < {prev_ds_cnt}")

        manifest = {
            "release": release_dir.name,
            "created_at": now_iso(),
            "mode": "incremental_parquet_sync",
            "datasets_count": int(new_ds_cnt),
            "publications_count": int(new_pub_cnt),
            "samples_sources": len(existing_sample_sources | touched_sample_sources),
            "artifacts_sources": len(existing_art_sources | touched_art_sources),
            "delta_datasets_ingested": int(delta_ds_cnt),
            "l1_deduped": int(l1_deduped),
            "rewritten_sample_partitions": sorted(touched_sample_sources),
            "rewritten_artifact_partitions": sorted(touched_art_sources),
        }
        (release_dir / "_manifest.json").write_text(json.dumps(manifest, indent=2))
        atomic_symlink_swap(release_dir, CURRENT_LINK)
        cleanup_old_releases(keep=max(1, keep_releases))
        log.info("Published atomic incremental release -> %s (%s)", CURRENT_LINK, json.dumps(manifest))
        return {"status": "published", **manifest}
    finally:
        con.close()
        shutil.rmtree(spill_dir, ignore_errors=True)


def main() -> int:
    p = argparse.ArgumentParser(description="Parquet-native tiered incremental sync & DuckDB dedup merger")
    p.add_argument("--status", action="store_true", help="Print schedule status for all configured sources")
    p.add_argument("--dry-run", action="store_true", help="Show which sources are due without fetching or writing")
    p.add_argument("--source", action="append", help="Restrict run to specific source(s), e.g. --source geo")
    p.add_argument("--force", action="store_true", help="Ignore interval_hours and force sync")
    p.add_argument("--max-items", type=int, default=None, help="Cap fetched items per source (useful for smoke tests)")
    p.add_argument("--delta-sqlite", type=str, default=None, help="Merge an external delta SQLite directly into Parquet")
    args = p.parse_args()

    policy = load_policy()
    state = load_state()
    defaults = policy.get("defaults") or {}
    memory_limit = str(defaults.get("memory_limit") or "2GB")
    keep_releases = int(defaults.get("keep_releases") or 2)
    overlap_days = int(defaults.get("overlap_days") or 1)

    if args.delta_sqlite:
        delta_path = Path(args.delta_sqlite).resolve()
        if not delta_path.exists():
            log.error("Delta SQLite not found: %s", delta_path)
            return 1
        res = merge_delta_into_parquet(delta_path, memory_limit=memory_limit, keep_releases=keep_releases)
        if args.source:
            for src in args.source:
                state.setdefault("sources", {})[src] = {"last_synced_at": now_iso(), "last_result": res}
            save_state(state)
        log.info("Delta merge result: %s", json.dumps(res, ensure_ascii=False))
        return 0

    all_sources = policy.get("sources") or {}
    target_sources = args.source if args.source else list(all_sources.keys())

    due_list: list[tuple[str, dict[str, Any], datetime]] = []
    for src in target_sources:
        cfg = all_sources.get(src)
        if not cfg:
            log.warning("Unknown source in policy: %s", src)
            continue
        src_state = state.get("sources", {}).get(src) or {}
        due, reason = is_due(cfg, src_state, force=args.force)
        last_at = src_state.get("last_synced_at") or state.get("bootstrapped_at") or now_iso()
        if args.status or args.dry_run:
            print(
                f"{src:16s} | tier={cfg.get('tier')} | every={cfg.get('interval_hours'):>4}h | "
                f"enabled={str(cfg.get('enabled')):5s} | due={str(due):5s} ({reason}) | last={last_at}"
            )
        if due and src in FETCHERS:
            since = (parse_iso(last_at) or now_dt()) - timedelta(days=overlap_days)
            due_list.append((src, cfg, since))

    if args.status or args.dry_run:
        return 0

    if not due_list:
        log.info("No sources due for incremental sync.")
        return 0

    tmp_delta = Path(tempfile.mkdtemp(prefix="omics_delta_")) / "delta.db"
    try:
        conn = init_delta_db(tmp_delta)
        synced_counts: dict[str, int] = {}
        for src, _cfg, since_dt in due_list:
            log.info("Fetching incremental delta for %s (since=%s) ...", src, since_dt.strftime("%Y-%m-%d"))
            cnt = FETCHERS[src](conn, since_dt, max_items=args.max_items)
            synced_counts[src] = cnt
            log.info("Fetched %s: %d datasets", src, cnt)
        conn.close()

        total_fetched = sum(synced_counts.values())
        if total_fetched == 0:
            log.info("All due sources returned 0 new datasets; advancing watermarks.")
            ts = now_iso()
            for src in synced_counts:
                state.setdefault("sources", {})[src] = {"last_synced_at": ts, "last_fetched": 0}
            save_state(state)
            return 0

        res = merge_delta_into_parquet(tmp_delta, memory_limit=memory_limit, keep_releases=keep_releases)
        ts = now_iso()
        for src, cnt in synced_counts.items():
            state.setdefault("sources", {})[src] = {
                "last_synced_at": ts,
                "last_fetched": cnt,
                "release": res.get("release"),
            }
        save_state(state)
        return 0
    finally:
        shutil.rmtree(tmp_delta.parent, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
