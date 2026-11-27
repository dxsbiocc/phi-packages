#!/usr/bin/env python3
"""Generic, idempotent merger from any omics_*_staging.db into omics_v1.db.

Features:
- Works across all staging DBs conforming to schema_v1.sql (datasets, samples,
  dataset_artifacts, dataset_members, publications, dataset_publications).
- L1 deduplication: resolves existing canonical_id in main via dataset_members
  (source, accession) or accession match, and resolves publications by DOI.
- Incremental watermark per staging DB stored in main.sync_cursors under
  source = 'merge:<db_name>'.
- Short-transaction batching with lock retry so main DB stays responsive.

Usage:
  python3 merge_staging_to_main.py --dry-run
  python3 merge_staging_to_main.py --db omics_pride_staging.db --batch-size 200
  python3 merge_staging_to_main.py --all-finished
"""

from __future__ import annotations

import argparse
import json
import logging
import sqlite3
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

try:
    from .paths import CODE_ROOT, data_dir, logs_dir
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from paths import CODE_ROOT, data_dir, logs_dir  # type: ignore

ROOT = CODE_ROOT
DATA_DIR = data_dir()
LOGS_DIR = logs_dir()
LOGS_DIR.mkdir(parents=True, exist_ok=True)
MAIN_DB = DATA_DIR / "omics_v1.db"
TZ = timezone(timedelta(hours=8))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOGS_DIR / "merge_staging_to_main.log"),
    ],
)
log = logging.getLogger("merge_staging_to_main")


def now_iso() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")


def connect_rw(path: Path, busy_ms: int = 60000) -> sqlite3.Connection:
    conn = sqlite3.connect(path, timeout=max(busy_ms / 1000, 10))
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute(f"PRAGMA busy_timeout = {busy_ms}")
    return conn


def connect_ro(path: Path) -> sqlite3.Connection:
    shm = path.with_name(path.name + "-shm")
    uri = f"file:{path}?mode=ro" if shm.exists() else f"file:{path}?mode=ro&immutable=1"
    conn = sqlite3.connect(uri, uri=True, timeout=15)
    conn.execute("PRAGMA temp_store = MEMORY")
    return conn


def is_staging_finished(st: sqlite3.Connection) -> tuple[bool, list[str]]:
    rows = st.execute("SELECT source, cursor_value FROM sync_cursors").fetchall()
    if not rows:
        return False, []
    sources = []
    all_done = True
    for src, val in rows:
        sources.append(src)
        try:
            obj = json.loads(val) if val else {}
        except Exception:
            obj = {}
        if not (bool(obj.get("done")) or bool(obj.get("completed_at"))):
            all_done = False
    return all_done, sources


def load_merge_cursor(main: sqlite3.Connection, db_name: str) -> dict[str, Any]:
    row = main.execute(
        "SELECT cursor_value FROM sync_cursors WHERE source = ?",
        (f"merge:{db_name}",),
    ).fetchone()
    if not row or not row[0]:
        return {"last_updated_at": "", "last_canonical_id": "", "merged_datasets": 0, "l1_deduped": 0}
    data = json.loads(row[0])
    data.setdefault("last_canonical_id", "")
    return data


def save_merge_cursor(
    main: sqlite3.Connection, db_name: str, cur: dict[str, Any], error: str | None = None
) -> None:
    main.execute(
        """
        INSERT INTO sync_cursors(source, cursor_value, cursor_kind, last_run_at, last_error)
        VALUES (?, ?, 'merge_watermark', ?, ?)
        ON CONFLICT(source) DO UPDATE SET
            cursor_value = excluded.cursor_value,
            cursor_kind  = excluded.cursor_kind,
            last_run_at  = excluded.last_run_at,
            last_error   = excluded.last_error
        """,
        (f"merge:{db_name}", json.dumps(cur), now_iso(), error),
    )


def resolve_canonical_id(
    main: sqlite3.Connection,
    st: sqlite3.Connection,
    staging_ds_id: str,
) -> tuple[str, bool]:
    """L1 dedup: check if any member (source, accession) or accession already exists in main."""
    members = st.execute(
        "SELECT source, accession FROM dataset_members WHERE dataset_id = ?",
        (staging_ds_id,),
    ).fetchall()
    for src, acc in members:
        row = main.execute(
            "SELECT m.dataset_id FROM dataset_members m JOIN datasets d ON d.canonical_id = m.dataset_id WHERE m.source = ? AND m.accession = ?",
            (src, acc),
        ).fetchone()
        if row:
            return row[0], row[0] != staging_ds_id
    for _src, acc in members:
        row = main.execute(
            "SELECT m.dataset_id FROM dataset_members m JOIN datasets d ON d.canonical_id = m.dataset_id WHERE m.accession = ? LIMIT 1",
            (acc,),
        ).fetchone()
        if row:
            return row[0], row[0] != staging_ds_id
    return staging_ds_id, False


def merge_dataset_batch(
    main: sqlite3.Connection,
    st: sqlite3.Connection,
    ds_rows: list[tuple],
    stats: dict[str, int],
) -> str:
    """Merge a batch of datasets and their child records inside one transaction."""
    max_ts = ""
    id_map: dict[str, str] = {}

    for row in ds_rows:
        (
            st_id, title, title_norm, organism, tissue, disease_tags,
            assay, platforms, n_samples, access, license_val, confidence,
            content_hash, first_seen_at, updated_at, probe_status, last_checked,
        ) = row
        if updated_at and updated_at > max_ts:
            max_ts = updated_at

        target_id, deduped = resolve_canonical_id(main, st, st_id)
        id_map[st_id] = target_id
        if deduped:
            stats["l1_deduped"] += 1

        if target_id == st_id:
            main.execute(
                """
                INSERT INTO datasets(
                    canonical_id, title, title_norm, organism, tissue, disease_tags,
                    assay, platforms, n_samples, access, license, confidence,
                    content_hash, first_seen_at, updated_at, probe_status, last_checked
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(canonical_id) DO UPDATE SET
                    title         = COALESCE(excluded.title, datasets.title),
                    title_norm    = COALESCE(excluded.title_norm, datasets.title_norm),
                    organism      = COALESCE(excluded.organism, datasets.organism),
                    tissue        = COALESCE(excluded.tissue, datasets.tissue),
                    disease_tags  = COALESCE(excluded.disease_tags, datasets.disease_tags),
                    assay         = COALESCE(excluded.assay, datasets.assay),
                    platforms     = COALESCE(excluded.platforms, datasets.platforms),
                    n_samples     = COALESCE(excluded.n_samples, datasets.n_samples),
                    access        = COALESCE(excluded.access, datasets.access),
                    license       = COALESCE(excluded.license, datasets.license),
                    confidence    = COALESCE(excluded.confidence, datasets.confidence),
                    content_hash  = excluded.content_hash,
                    updated_at    = excluded.updated_at,
                    probe_status  = excluded.probe_status,
                    last_checked  = excluded.last_checked
                """,
                row,
            )
        else:
            # Enrich existing canonical dataset without overwriting primary fields
            main.execute(
                """
                UPDATE datasets SET
                    organism     = COALESCE(organism, ?),
                    tissue       = COALESCE(tissue, ?),
                    disease_tags = COALESCE(disease_tags, ?),
                    assay        = COALESCE(assay, ?),
                    platforms    = COALESCE(platforms, ?),
                    n_samples    = COALESCE(n_samples, ?),
                    updated_at   = COALESCE(?, updated_at),
                    last_checked = COALESCE(?, last_checked)
                WHERE canonical_id = ?
                """,
                (organism, tissue, disease_tags, assay, platforms, n_samples, updated_at, last_checked, target_id),
            )
        stats["datasets"] += 1

    st_ids = list(id_map.keys())
    placeholders = ",".join("?" * len(st_ids))

    # 1. dataset_members
    for ds_id, src, acc, role, ext_url in st.execute(
        f"SELECT dataset_id, source, accession, role, external_url FROM dataset_members WHERE dataset_id IN ({placeholders})",
        st_ids,
    ):
        main.execute(
            """
            INSERT INTO dataset_members(dataset_id, source, accession, role, external_url)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(source, accession) DO NOTHING
            """,
            (id_map[ds_id], src, acc, role, ext_url),
        )
        stats["members"] += 1

    # 2. samples
    for srow in st.execute(
        f"""
        SELECT sample_id, dataset_id, title, organism, tissue, disease, condition,
               sex, age, attributes, run_ids, content_hash, updated_at
        FROM samples WHERE dataset_id IN ({placeholders})
        """,
        st_ids,
    ):
        s_list = list(srow)
        s_list[1] = id_map[s_list[1]]
        main.execute(
            """
            INSERT INTO samples(
                sample_id, dataset_id, title, organism, tissue, disease, condition,
                sex, age, attributes, run_ids, content_hash, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(sample_id) DO UPDATE SET
                dataset_id   = excluded.dataset_id,
                title        = COALESCE(excluded.title, samples.title),
                organism     = COALESCE(excluded.organism, samples.organism),
                tissue       = COALESCE(excluded.tissue, samples.tissue),
                disease      = COALESCE(excluded.disease, samples.disease),
                condition    = COALESCE(excluded.condition, samples.condition),
                sex          = COALESCE(excluded.sex, samples.sex),
                age          = COALESCE(excluded.age, samples.age),
                attributes   = COALESCE(excluded.attributes, samples.attributes),
                run_ids      = COALESCE(excluded.run_ids, samples.run_ids),
                content_hash = excluded.content_hash,
                updated_at   = excluded.updated_at
            """,
            s_list,
        )
        stats["samples"] += 1

    # 3. dataset_artifacts
    for arow in st.execute(
        f"""
        SELECT artifact_id, dataset_id, sample_id, kind, role, title, filename,
               format, uri, size_bytes, access, license, source, derived_from,
               probe_status, last_checked, content_hash
        FROM dataset_artifacts WHERE dataset_id IN ({placeholders})
        """,
        st_ids,
    ):
        a_list = list(arow)
        a_list[1] = id_map[a_list[1]]
        if a_list[2] is not None:
            has_sample = main.execute("SELECT 1 FROM samples WHERE sample_id = ?", (a_list[2],)).fetchone()
            if not has_sample:
                a_list[2] = None
        main.execute(
            """
            INSERT INTO dataset_artifacts(
                artifact_id, dataset_id, sample_id, kind, role, title, filename,
                format, uri, size_bytes, access, license, source, derived_from,
                probe_status, last_checked, content_hash
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(artifact_id) DO UPDATE SET
                dataset_id   = excluded.dataset_id,
                uri          = excluded.uri,
                size_bytes   = COALESCE(excluded.size_bytes, dataset_artifacts.size_bytes),
                last_checked = excluded.last_checked,
                content_hash = excluded.content_hash
            """,
            a_list,
        )
        stats["artifacts"] += 1

    # 4. publications & dataset_publications (with DOI conflict resolution)
    pmid_map: dict[str, str] = {}
    for prow in st.execute(
        f"""
        SELECT DISTINCT p.pmid, p.doi, p.title, p.authors, p.journal, p.year, p.abstract, p.url, p.fetched_at
        FROM publications p
        JOIN dataset_publications dp ON dp.pmid = p.pmid
        WHERE dp.dataset_id IN ({placeholders})
        """,
        st_ids,
    ):
        pmid, doi = prow[0], prow[1]
        target_pmid = pmid
        if doi:
            existing_doi = main.execute(
                "SELECT pmid FROM publications WHERE doi = ?", (doi,)
            ).fetchone()
            if existing_doi:
                target_pmid = existing_doi[0]
        pmid_map[pmid] = target_pmid
        if target_pmid == pmid:
            main.execute(
                """
                INSERT INTO publications(pmid, doi, title, authors, journal, year, abstract, url, fetched_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(pmid) DO UPDATE SET
                    doi        = COALESCE(excluded.doi, publications.doi),
                    title      = COALESCE(excluded.title, publications.title),
                    authors    = COALESCE(excluded.authors, publications.authors),
                    journal    = COALESCE(excluded.journal, publications.journal),
                    year       = COALESCE(excluded.year, publications.year),
                    abstract   = COALESCE(excluded.abstract, publications.abstract),
                    url        = COALESCE(excluded.url, publications.url),
                    fetched_at = excluded.fetched_at
                """,
                prow,
            )
            stats["publications"] += 1

    for ds_id, pmid, relation, source in st.execute(
        f"SELECT dataset_id, pmid, relation, source FROM dataset_publications WHERE dataset_id IN ({placeholders})",
        st_ids,
    ):
        main.execute(
            """
            INSERT INTO dataset_publications(dataset_id, pmid, relation, source)
            VALUES (?, ?, ?, ?)
            ON CONFLICT(dataset_id, pmid) DO NOTHING
            """,
            (id_map[ds_id], pmid_map.get(pmid, pmid), relation, source),
        )

    return max_ts


def merge_staging_db(
    staging_path: Path,
    batch_size: int,
    dry_run: bool,
    include_unfinished: bool,
) -> dict[str, Any]:
    st = connect_ro(staging_path)
    finished, sources = is_staging_finished(st)
    if not finished and not include_unfinished:
        st.close()
        return {"db": staging_path.name, "status": "skipped_unfinished", "sources": sources}

    main_ro = connect_ro(MAIN_DB)
    cur = load_merge_cursor(main_ro, staging_path.name)
    main_ro.close()

    watermark = cur.get("last_updated_at") or ""
    last_cid = cur.get("last_canonical_id") or ""
    pending_count = st.execute(
        "SELECT COUNT(*) FROM datasets WHERE COALESCE(updated_at, '') > ? OR (COALESCE(updated_at, '') = ? AND canonical_id > ?)",
        (watermark, watermark, last_cid),
    ).fetchone()[0]

    if dry_run or pending_count == 0:
        st.close()
        return {
            "db": staging_path.name,
            "status": "dry_run" if dry_run else "up_to_date",
            "finished": finished,
            "watermark": watermark,
            "pending_datasets": pending_count,
        }

    stats = {"datasets": 0, "samples": 0, "artifacts": 0, "members": 0, "publications": 0, "l1_deduped": 0}
    last_ts = watermark

    while True:
        rows = st.execute(
            """
            SELECT canonical_id, title, title_norm, organism, tissue, disease_tags,
                   assay, platforms, n_samples, access, license, confidence,
                   content_hash, first_seen_at, updated_at, probe_status, last_checked
            FROM datasets
            WHERE COALESCE(updated_at, '') > ?
               OR (COALESCE(updated_at, '') = ? AND canonical_id > ?)
            ORDER BY COALESCE(updated_at, '') ASC, canonical_id ASC
            LIMIT ?
            """,
            (last_ts, last_ts, last_cid, batch_size),
        ).fetchall()
        if not rows:
            break

        for attempt in range(1, 20):
            main = connect_rw(MAIN_DB)
            try:
                merge_dataset_batch(main, st, rows, stats)
                last_ts = rows[-1][14] or ""
                last_cid = rows[-1][0]
                cur["last_updated_at"] = last_ts
                cur["last_canonical_id"] = last_cid
                cur["merged_datasets"] = cur.get("merged_datasets", 0) + len(rows)
                cur["l1_deduped"] = cur.get("l1_deduped", 0) + stats["l1_deduped"]
                save_merge_cursor(main, staging_path.name, cur)
                main.commit()
                break
            except sqlite3.OperationalError as e:
                main.rollback()
                if "locked" in str(e).lower() and attempt < 19:
                    time.sleep(min(2 * attempt, 20))
                    continue
                raise
            finally:
                main.close()

        log.info("%s progress: merged=%d/%d l1_deduped=%d", staging_path.name, stats["datasets"], pending_count, stats["l1_deduped"])
        if len(rows) < batch_size:
            break

    st.close()
    return {"db": staging_path.name, "status": "merged", "watermark": last_ts, **stats}


def main() -> int:
    p = argparse.ArgumentParser(description="Generic staging -> omics_v1.db merger")
    p.add_argument("--db", action="append", help="Specific staging DB filename(s) to merge")
    p.add_argument("--all-finished", action="store_true", help="Merge all staging DBs marked done in sync_cursors")
    p.add_argument("--include-unfinished", action="store_true", help="Also merge staging DBs not yet marked done")
    p.add_argument("--batch-size", type=int, default=200, help="Datasets per transaction (default: 200)")
    p.add_argument("--dry-run", action="store_true", help="Report pending counts without writing omics_v1.db")
    args = p.parse_args()

    if args.db:
        targets = [DATA_DIR / name if not Path(name).is_absolute() else Path(name) for name in args.db]
    else:
        targets = sorted(DATA_DIR.glob("omics_*_staging.db"))

    m = connect_rw(MAIN_DB)
    m.execute("CREATE INDEX IF NOT EXISTS idx_members_accession ON dataset_members(accession)")
    m.commit()
    m.close()

    if not targets:
        log.warning("No staging databases found.")
        return 0

    results = []
    for path in targets:
        if not path.exists():
            log.warning("Missing staging DB: %s", path)
            continue
        res = merge_staging_db(
            path,
            batch_size=max(1, args.batch_size),
            dry_run=args.dry_run,
            include_unfinished=args.include_unfinished or bool(args.db),
        )
        results.append(res)
        log.info("Result: %s", json.dumps(res, ensure_ascii=False))

    return 0


if __name__ == "__main__":
    sys.exit(main())
