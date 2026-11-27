#!/usr/bin/env python3
"""
CELLxGENE Discover / Curation catalog ingest (scan / staging phase).

Writes ONLY to omics_cellxgene_staging.db — never to omics_v1.db.
Catalog-only: dataset cards + asset URL pointers (no h5ad download).
Resumable via sync_cursors (source=cellxgene_full_scan).
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
STAGING_DB = _DATA_DIR / "omics_cellxgene_staging.db"
SCHEMA = SCHEMA_PATH
API = "https://api.cellxgene.cziscience.com/curation/v1"
CURSOR_SOURCE = "cellxgene_full_scan"
MIN_INTERVAL = 0.30
MAX_RETRIES = 8
TZ = timezone(timedelta(hours=8))
SAVE_EVERY = 50
LOG_EVERY = 100
BATCH_COMMIT = 25
# Cap cell_type labels stored per dataset to keep DB lean
MAX_CELL_TYPES = 80

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_LOGS_DIR / "cellxgene_full_scan.log"),
    ],
)
log = logging.getLogger("cellxgene_full_scan")
_last_req = 0.0


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def throttle() -> None:
    global _last_req
    elapsed = time.monotonic() - _last_req
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_req = time.monotonic()


def cxg_get(url: str, params: dict | None = None, timeout: Any = TIMEOUT) -> Any:
    last_err: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        throttle()
        try:
            r = http_get(url, params=params, timeout=timeout)
            if r.status_code == 429 or r.status_code >= 500:
                wait = min(2**attempt, 30)
                log.warning("HTTP %s %s (attempt %s), sleep %ss", r.status_code, url, attempt, wait)
                time.sleep(wait)
                continue
            if r.status_code in (400, 404):
                return None
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last_err = e
            wait = min(2**attempt, 30)
            log.warning("Request failed attempt %s %s: %s; sleep %ss", attempt, url, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"CELLxGENE request failed after {MAX_RETRIES}: {last_err}")


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
        CREATE TABLE IF NOT EXISTS src_cellxgene_dataset (
            dataset_id             TEXT PRIMARY KEY,  -- UUID
            dataset_version_id     TEXT,
            title                  TEXT,
            collection_id          TEXT,
            collection_name        TEXT,
            collection_doi         TEXT,
            organisms              TEXT,   -- JSON labels
            tissues                TEXT,
            diseases               TEXT,
            assays                 TEXT,
            cell_types             TEXT,   -- capped JSON labels
            cell_count             INTEGER,
            primary_cell_count     INTEGER,
            donor_count            INTEGER,
            mean_genes_per_cell    REAL,
            is_primary_data        TEXT,
            schema_version         TEXT,
            published_at           TEXT,
            revised_at             TEXT,
            citation               TEXT,
            explorer_url           TEXT,
            assets                 TEXT,   -- JSON [{filetype,filesize,url}]
            spatial                TEXT,
            suspension_type        TEXT,
            sex                    TEXT,
            development_stage      TEXT,
            visibility             TEXT,
            raw_json               TEXT,   -- slim snapshot (no full cell_type dump)
            fetched_at             TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_src_cxg_fetched_at ON src_cellxgene_dataset(fetched_at);
        CREATE INDEX IF NOT EXISTS idx_src_cxg_collection ON src_cellxgene_dataset(collection_id);
        CREATE TABLE IF NOT EXISTS cellxgene_merge_state (
            dataset_id TEXT PRIMARY KEY,
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
            "skipped_tombstone": 0,
            "ingested_new": 0,
            "src_upserted": 0,
            "artifacts_created": 0,
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
        (CURSOR_SOURCE, json.dumps(cur), "index", now_iso(), error),
    )
    conn.commit()


def labels(items: Any, limit: int | None = None) -> list[str]:
    out: list[str] = []
    if not items:
        return out
    if isinstance(items, list):
        for it in items:
            if isinstance(it, dict):
                lab = it.get("label") or it.get("name")
                if lab:
                    out.append(str(lab))
            elif it:
                out.append(str(it))
    elif isinstance(items, str):
        out.append(items)
    out = list(dict.fromkeys(out))
    if limit is not None:
        out = out[:limit]
    return out


def join_labels(items: Any, sep: str = "; ", limit: int | None = None) -> str | None:
    labs = labels(items, limit=limit)
    return sep.join(labs) if labs else None


def as_int(val: Any) -> int | None:
    if val is None or val == "":
        return None
    try:
        return int(val)
    except (TypeError, ValueError):
        return None


def as_float(val: Any) -> float | None:
    if val is None or val == "":
        return None
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def slim_dataset(d: dict) -> dict:
    """Keep catalog fields; drop bulky nested ontology dumps beyond labels."""
    assets = []
    for a in d.get("assets") or []:
        if not isinstance(a, dict):
            continue
        assets.append(
            {
                "filetype": a.get("filetype"),
                "filesize": a.get("filesize"),
                "url": a.get("url"),
            }
        )
    return {
        "dataset_id": d.get("dataset_id"),
        "dataset_version_id": d.get("dataset_version_id"),
        "title": d.get("title"),
        "collection_id": d.get("collection_id"),
        "collection_name": d.get("collection_name"),
        "collection_doi": d.get("collection_doi"),
        "collection_doi_label": d.get("collection_doi_label"),
        "organism": labels(d.get("organism")),
        "tissue": labels(d.get("tissue")),
        "disease": labels(d.get("disease")),
        "assay": labels(d.get("assay")),
        "cell_type": labels(d.get("cell_type"), limit=MAX_CELL_TYPES),
        "cell_count": d.get("cell_count"),
        "primary_cell_count": d.get("primary_cell_count"),
        "donor_id_n": len(d.get("donor_id") or []) if isinstance(d.get("donor_id"), list) else None,
        "mean_genes_per_cell": d.get("mean_genes_per_cell"),
        "is_primary_data": d.get("is_primary_data"),
        "schema_version": d.get("schema_version"),
        "published_at": d.get("published_at"),
        "revised_at": d.get("revised_at"),
        "citation": (d.get("citation") or "")[:2000] or None,
        "explorer_url": d.get("explorer_url"),
        "assets": assets,
        "spatial": d.get("spatial"),
        "suspension_type": d.get("suspension_type"),
        "sex": labels(d.get("sex")),
        "development_stage": labels(d.get("development_stage"), limit=40),
        "visibility": d.get("visibility"),
        "tombstone": d.get("tombstone"),
        "processing_status": d.get("processing_status"),
    }


def fetch_catalog() -> list[dict]:
    log.info("Fetching CELLxGENE curation datasets catalog...")
    data = cxg_get(f"{API}/datasets", timeout=(15, 180))
    if not isinstance(data, list):
        raise RuntimeError(f"Unexpected CELLxGENE datasets payload: {type(data)}")
    log.info("Catalog has %s datasets", len(data))
    return data


def staging_ingest(conn: sqlite3.Connection, d: dict, cur: dict) -> None:
    dataset_id = (d.get("dataset_id") or "").strip()
    if not dataset_id:
        cur["errors"] = cur.get("errors", 0) + 1
        log.warning("skip dataset without id")
        return
    if d.get("tombstone"):
        cur["skipped_tombstone"] = cur.get("skipped_tombstone", 0) + 1
        return

    exists = conn.execute(
        "SELECT 1 FROM src_cellxgene_dataset WHERE dataset_id=? LIMIT 1", (dataset_id,)
    ).fetchone()
    if exists:
        cur["skipped_existing"] = cur.get("skipped_existing", 0) + 1
        return

    fetched_at = now_iso()
    slim = slim_dataset(d)
    title = slim.get("title") or dataset_id
    organisms = slim.get("organism") or []
    tissues = slim.get("tissue") or []
    diseases = slim.get("disease") or []
    assays = slim.get("assay") or []
    cell_types = slim.get("cell_type") or []
    assets = slim.get("assets") or []
    cell_count = as_int(slim.get("cell_count"))
    primary_cell_count = as_int(slim.get("primary_cell_count"))
    donor_count = slim.get("donor_id_n")
    is_primary = slim.get("is_primary_data")
    if isinstance(is_primary, list):
        is_primary = json.dumps(is_primary)
    elif is_primary is not None:
        is_primary = str(is_primary)

    conn.execute(
        """
        INSERT INTO src_cellxgene_dataset(
            dataset_id, dataset_version_id, title, collection_id, collection_name,
            collection_doi, organisms, tissues, diseases, assays, cell_types,
            cell_count, primary_cell_count, donor_count, mean_genes_per_cell,
            is_primary_data, schema_version, published_at, revised_at, citation,
            explorer_url, assets, spatial, suspension_type, sex, development_stage,
            visibility, raw_json, fetched_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(dataset_id) DO UPDATE SET
            dataset_version_id=excluded.dataset_version_id, title=excluded.title,
            collection_id=excluded.collection_id, collection_name=excluded.collection_name,
            collection_doi=excluded.collection_doi, organisms=excluded.organisms,
            tissues=excluded.tissues, diseases=excluded.diseases, assays=excluded.assays,
            cell_types=excluded.cell_types, cell_count=excluded.cell_count,
            primary_cell_count=excluded.primary_cell_count, donor_count=excluded.donor_count,
            mean_genes_per_cell=excluded.mean_genes_per_cell,
            is_primary_data=excluded.is_primary_data, schema_version=excluded.schema_version,
            published_at=excluded.published_at, revised_at=excluded.revised_at,
            citation=excluded.citation, explorer_url=excluded.explorer_url,
            assets=excluded.assets, spatial=excluded.spatial,
            suspension_type=excluded.suspension_type, sex=excluded.sex,
            development_stage=excluded.development_stage, visibility=excluded.visibility,
            raw_json=excluded.raw_json, fetched_at=excluded.fetched_at
        """,
        (
            dataset_id,
            slim.get("dataset_version_id"),
            title,
            slim.get("collection_id"),
            slim.get("collection_name"),
            slim.get("collection_doi"),
            json.dumps(organisms, ensure_ascii=False),
            json.dumps(tissues, ensure_ascii=False),
            json.dumps(diseases, ensure_ascii=False),
            json.dumps(assays, ensure_ascii=False),
            json.dumps(cell_types, ensure_ascii=False),
            cell_count,
            primary_cell_count,
            donor_count,
            as_float(slim.get("mean_genes_per_cell")),
            is_primary,
            slim.get("schema_version"),
            slim.get("published_at"),
            slim.get("revised_at"),
            slim.get("citation"),
            slim.get("explorer_url"),
            json.dumps(assets, ensure_ascii=False),
            json.dumps(slim.get("spatial"), ensure_ascii=False) if slim.get("spatial") is not None else None,
            json.dumps(slim.get("suspension_type"), ensure_ascii=False)
            if slim.get("suspension_type") is not None
            else None,
            json.dumps(slim.get("sex") or [], ensure_ascii=False),
            json.dumps(slim.get("development_stage") or [], ensure_ascii=False),
            slim.get("visibility"),
            json.dumps(slim, ensure_ascii=False, default=str),
            fetched_at,
        ),
    )
    cur["src_upserted"] = cur.get("src_upserted", 0) + 1

    canonical_id = f"cellxgene:{dataset_id}"
    ch = content_hash(
        {
            "id": dataset_id,
            "ver": slim.get("dataset_version_id"),
            "title": title,
            "cells": cell_count,
            "n_assets": len(assets),
        }
    )
    existing = conn.execute(
        "SELECT first_seen_at FROM datasets WHERE canonical_id=?", (canonical_id,)
    ).fetchone()
    first_seen = existing[0] if existing else fetched_at
    organism = "; ".join(organisms) if organisms else None
    tissue = "; ".join(tissues) if tissues else None
    disease_tags = json.dumps(diseases, ensure_ascii=False) if diseases else None
    assay = "; ".join(assays) if assays else "single-cell"
    platforms = json.dumps(assays, ensure_ascii=False) if assays else None

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
            cell_count,
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

    collection_url = None
    if slim.get("collection_id"):
        collection_url = f"https://cellxgene.cziscience.com/collections/{slim['collection_id']}"
    explorer = slim.get("explorer_url") or (
        f"https://cellxgene.cziscience.com/e/{slim.get('dataset_version_id')}.cxg/"
        if slim.get("dataset_version_id")
        else None
    )
    conn.execute(
        """
        INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
        VALUES (?, 'cellxgene', ?, 'dataset', ?)
        ON CONFLICT(source, accession) DO UPDATE SET
            dataset_id=excluded.dataset_id, role=excluded.role, external_url=excluded.external_url
        """,
        (canonical_id, dataset_id, explorer or collection_url),
    )
    if slim.get("collection_id"):
        conn.execute(
            """
            INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
            VALUES (?, 'cellxgene_collection', ?, 'collection', ?)
            ON CONFLICT(source, accession) DO NOTHING
            """,
            (canonical_id, slim["collection_id"], collection_url),
        )

    # Optional DOI publication link
    doi = (slim.get("collection_doi") or "").strip()
    if doi:
        # strip URL prefix if present
        doi = re.sub(r"^https?://(dx\.)?doi\.org/", "", doi, flags=re.I)
        pk = f"DOI:{doi}"
        try:
            conn.execute(
                """
                INSERT INTO publications(pmid, doi, title, authors, journal, year, abstract, url, fetched_at)
                VALUES (?, ?, ?, NULL, NULL, NULL, NULL, ?, ?)
                ON CONFLICT(pmid) DO UPDATE SET
                    doi=COALESCE(excluded.doi, publications.doi),
                    fetched_at=excluded.fetched_at
                """,
                (pk, doi, slim.get("collection_name"), f"https://doi.org/{doi}", fetched_at),
            )
            conn.execute(
                """
                INSERT INTO dataset_publications(dataset_id, pmid, relation, source)
                VALUES (?, ?, 'collection', 'cellxgene')
                ON CONFLICT(dataset_id, pmid) DO NOTHING
                """,
                (canonical_id, pk),
            )
        except sqlite3.IntegrityError:
            row = conn.execute(
                "SELECT pmid FROM publications WHERE doi=? LIMIT 1", (doi,)
            ).fetchone()
            if row:
                conn.execute(
                    """
                    INSERT INTO dataset_publications(dataset_id, pmid, relation, source)
                    VALUES (?, ?, 'collection', 'cellxgene')
                    ON CONFLICT(dataset_id, pmid) DO NOTHING
                    """,
                    (canonical_id, row[0]),
                )

    conn.execute("DELETE FROM dataset_artifacts WHERE dataset_id=?", (canonical_id,))
    for a in assets:
        url = a.get("url")
        if not url:
            continue
        ft = (a.get("filetype") or "FILE").upper()
        kind = "processed" if ft in ("H5AD", "RDS", "CXG") else "supplementary"
        role = ft.lower()
        fname = url.rstrip("/").rsplit("/", 1)[-1]
        artifact_id = f"cellxgene:{dataset_id}:{hashlib.sha1(url.encode()).hexdigest()[:16]}"
        ach = content_hash({"uri": url, "ft": ft, "sz": a.get("filesize")})
        size_i = as_int(a.get("filesize"))
        conn.execute(
            """
            INSERT INTO dataset_artifacts(
                artifact_id, dataset_id, sample_id, kind, role, title, filename,
                format, uri, size_bytes, access, license, source, derived_from,
                probe_status, last_checked, content_hash
            ) VALUES (?, ?, NULL, ?, ?, ?, ?, ?, ?, ?, 'public', NULL, 'cellxgene', NULL,
                      'listed', ?, ?)
            ON CONFLICT(artifact_id) DO UPDATE SET
                uri=excluded.uri, size_bytes=excluded.size_bytes,
                kind=excluded.kind, role=excluded.role,
                last_checked=excluded.last_checked, content_hash=excluded.content_hash
            """,
            (
                artifact_id,
                canonical_id,
                kind,
                role,
                f"{title} ({ft})",
                fname,
                ft.lower(),
                url,
                size_i,
                fetched_at,
                ach,
            ),
        )
        cur["artifacts_created"] = cur.get("artifacts_created", 0) + 1


def main() -> int:
    log.info("Starting CELLxGENE SCAN-ONLY (staging=%s)", STAGING_DB)
    init_staging()
    st = connect(STAGING_DB, busy_ms=5000)
    cur = load_cursor(st)

    try:
        catalog = fetch_catalog()
        total = len(catalog)
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
            d = catalog[idx]
            try:
                staging_ingest(st, d if isinstance(d, dict) else {}, cur)
                pending_commit += 1
            except Exception as e:
                cur["errors"] = cur.get("errors", 0) + 1
                did = (d or {}).get("dataset_id") if isinstance(d, dict) else None
                log.exception("Failed staging %s: %s", did, e)
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
                    "progress index=%s/%s scanned=%s new=%s skipped=%s src=%s artifacts=%s errors=%s",
                    idx,
                    total,
                    cur.get("scanned"),
                    cur.get("ingested_new"),
                    cur.get("skipped_existing"),
                    cur.get("src_upserted"),
                    cur.get("artifacts_created"),
                    cur.get("errors"),
                )

        if pending_commit:
            st.commit()
        cur["index"] = total
        cur["completed_at"] = now_iso()
        save_cursor(st, cur)
        log.info("Finished CELLxGENE scan: %s", json.dumps(cur))
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
