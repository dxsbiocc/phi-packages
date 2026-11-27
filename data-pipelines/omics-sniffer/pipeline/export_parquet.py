#!/usr/bin/env python3
"""Atomic, source-partitioned Parquet exporter for omics_v1.db -> serving/current.

Design principles:
- SQLite (omics_v1.db) remains the Single Source of Truth (SSOT) for writes.
- Parquet under serving/ is a read-only, immutable OLAP layer for Phi MCP.
- Excludes raw src_* mirror tables (XML/JSON blobs).
- Compact tables (datasets, publications, dataset_publications, dataset_members)
  are exported directly when changed.
- Large tables (samples, dataset_artifacts) are partitioned by `source` and reuse
  unchanged partitions via hardlink (zero-copy) based on per-source watermarks.
- Publishes atomically via staging directory + symlink swap (serving/current).

Requires: duckdb (pip install duckdb)
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import sqlite3
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

try:
    from .paths import CODE_ROOT, data_dir, serving_dir
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from paths import CODE_ROOT, data_dir, serving_dir  # type: ignore

ROOT = CODE_ROOT
MAIN_DB = data_dir() / "omics_v1.db"
SERVING_DIR = serving_dir()
RELEASES_DIR = SERVING_DIR / "releases"
CURRENT_LINK = SERVING_DIR / "current"
STATE_FILE = SERVING_DIR / ".export_state.json"
TZ = timezone(timedelta(hours=8))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("export_parquet")


def now_tag() -> str:
    return datetime.now(TZ).strftime("v_%Y%m%d_%H%M%S")


def load_state() -> dict[str, Any]:
    if STATE_FILE.exists():
        return json.loads(STATE_FILE.read_text())
    return {"tables": {}, "samples_by_source": {}, "artifacts_by_source": {}}


def save_state(state: dict[str, Any]) -> None:
    tmp = STATE_FILE.with_suffix(".tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True))
    os.replace(tmp, STATE_FILE)


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


def get_watermarks(db_path: Path) -> dict[str, Any]:
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True, timeout=15)
    wm: dict[str, Any] = {
        "datasets": conn.execute("SELECT COUNT(*), COALESCE(MAX(updated_at), '') FROM datasets").fetchone(),
        "publications": conn.execute("SELECT COUNT(*), COALESCE(MAX(fetched_at), '') FROM publications").fetchone(),
        "dataset_publications": conn.execute("SELECT COUNT(*) FROM dataset_publications").fetchone()[0],
        "dataset_members": conn.execute("SELECT COUNT(*) FROM dataset_members").fetchone()[0],
    }

    # Per-source watermarks for samples (derived from prefix before ':' or 'geo')
    samples_src = conn.execute(
        """
        SELECT
            CASE WHEN instr(sample_id, ':') > 0
                 THEN substr(sample_id, 1, instr(sample_id, ':') - 1)
                 ELSE 'geo'
            END AS src,
            COUNT(*) AS cnt,
            COALESCE(MAX(updated_at), '') AS max_ts
        FROM samples
        GROUP BY src
        """
    ).fetchall()
    wm["samples_by_source"] = {r[0]: [r[1], r[2]] for r in samples_src}

    # Per-source watermarks for dataset_artifacts
    art_src = conn.execute(
        """
        SELECT COALESCE(NULLIF(source, ''), 'unknown') AS src,
               COUNT(*) AS cnt,
               COALESCE(MAX(last_checked), '') AS max_ts
        FROM dataset_artifacts
        GROUP BY src
        """
    ).fetchall()
    wm["artifacts_by_source"] = {r[0]: [r[1], r[2]] for r in art_src}
    conn.close()
    return wm


def run_export(full: bool, dry_run: bool, keep_releases: int) -> int:
    try:
        import duckdb  # type: ignore
    except ImportError:
        log.error("Missing dependency: duckdb. Install via `pip install duckdb` or use a Python env with duckdb.")
        return 2

    if not MAIN_DB.exists():
        log.error("Main DB not found: %s", MAIN_DB)
        return 1

    SERVING_DIR.mkdir(parents=True, exist_ok=True)
    RELEASES_DIR.mkdir(parents=True, exist_ok=True)

    prev_state = {} if full else load_state()
    wm = get_watermarks(MAIN_DB)

    if dry_run:
        log.info("Dry-run watermarks: %s", json.dumps({
            "datasets": wm["datasets"],
            "publications": wm["publications"],
            "samples_sources": len(wm["samples_by_source"]),
            "artifacts_sources": len(wm["artifacts_by_source"]),
        }))
        return 0

    release_dir = RELEASES_DIR / now_tag()
    release_dir.mkdir(parents=True, exist_ok=True)
    prev_release = CURRENT_LINK.resolve() if CURRENT_LINK.exists() else None

    con = duckdb.connect(database=":memory:")
    con.execute("INSTALL sqlite; LOAD sqlite;")
    con.execute(f"ATTACH '{MAIN_DB}' AS omics (TYPE SQLITE, READ_ONLY);")

    # 1. Export compact tables (reuse if watermark unchanged and previous file exists)
    compact_queries = {
        "datasets": (
            wm["datasets"],
            "SELECT * FROM omics.datasets",
        ),
        "publications": (
            wm["publications"],
            "SELECT * FROM omics.publications",
        ),
        "dataset_publications": (
            wm["dataset_publications"],
            "SELECT * FROM omics.dataset_publications",
        ),
        "dataset_members": (
            wm["dataset_members"],
            "SELECT * FROM omics.dataset_members",
        ),
    }

    for tbl, (cur_sig, sql) in compact_queries.items():
        dst = release_dir / f"{tbl}.parquet"
        prev_sig = prev_state.get("tables", {}).get(tbl)
        prev_file = (prev_release / f"{tbl}.parquet") if prev_release else None
        if not full and prev_sig == list(cur_sig) if isinstance(cur_sig, tuple) else prev_sig == cur_sig:
            if prev_file and prev_file.exists():
                link_or_copy(prev_file, dst)
                log.info("Reused %s.parquet via hardlink", tbl)
                continue
        log.info("Exporting %s.parquet ...", tbl)
        con.execute(f"COPY ({sql}) TO '{dst}' (FORMAT PARQUET, COMPRESSION 'ZSTD', ROW_GROUP_SIZE 100000);")

    # 2. Export samples partitioned by source
    samples_dir = release_dir / "samples"
    samples_dir.mkdir(parents=True, exist_ok=True)
    for src, sig in wm["samples_by_source"].items():
        part_dir = samples_dir / f"source={src}"
        dst = part_dir / "data.parquet"
        prev_sig = prev_state.get("samples_by_source", {}).get(src)
        prev_file = (prev_release / "samples" / f"source={src}" / "data.parquet") if prev_release else None
        if not full and prev_sig == sig and prev_file and prev_file.exists():
            link_or_copy(prev_file, dst)
            log.info("Reused samples/source=%s", src)
            continue
        part_dir.mkdir(parents=True, exist_ok=True)
        log.info("Exporting samples/source=%s (rows=%d) ...", src, sig[0])
        cond = (
            "instr(sample_id, ':') = 0"
            if src == "geo"
            else f"split_part(sample_id, ':', 1) = '{src}'"
        )
        con.execute(
            f"COPY (SELECT * FROM omics.samples WHERE {cond}) "
            f"TO '{dst}' (FORMAT PARQUET, COMPRESSION 'ZSTD', ROW_GROUP_SIZE 100000);"
        )

    # 3. Export dataset_artifacts partitioned by source
    art_dir = release_dir / "dataset_artifacts"
    art_dir.mkdir(parents=True, exist_ok=True)
    for src, sig in wm["artifacts_by_source"].items():
        part_dir = art_dir / f"source={src}"
        dst = part_dir / "data.parquet"
        prev_sig = prev_state.get("artifacts_by_source", {}).get(src)
        prev_file = (prev_release / "dataset_artifacts" / f"source={src}" / "data.parquet") if prev_release else None
        if not full and prev_sig == sig and prev_file and prev_file.exists():
            link_or_copy(prev_file, dst)
            log.info("Reused dataset_artifacts/source=%s", src)
            continue
        part_dir.mkdir(parents=True, exist_ok=True)
        log.info("Exporting dataset_artifacts/source=%s (rows=%d) ...", src, sig[0])
        cond = (
            "(source IS NULL OR source = '')"
            if src == "unknown"
            else f"source = '{src}'"
        )
        con.execute(
            f"COPY (SELECT * FROM omics.dataset_artifacts WHERE {cond}) "
            f"TO '{dst}' (FORMAT PARQUET, COMPRESSION 'ZSTD', ROW_GROUP_SIZE 100000);"
        )

    # 4. Validate row counts via DuckDB before atomic swap
    ds_cnt = con.execute(f"SELECT COUNT(*) FROM read_parquet('{release_dir / 'datasets.parquet'}')").fetchone()[0]
    if ds_cnt != wm["datasets"][0]:
        raise RuntimeError(f"Validation failed for datasets: parquet={ds_cnt} != sqlite={wm['datasets'][0]}")

    manifest = {
        "release": release_dir.name,
        "created_at": datetime.now(TZ).isoformat(),
        "datasets_count": ds_cnt,
        "publications_count": wm["publications"][0],
        "samples_sources": len(wm["samples_by_source"]),
        "artifacts_sources": len(wm["artifacts_by_source"]),
    }
    (release_dir / "_manifest.json").write_text(json.dumps(manifest, indent=2))
    con.close()

    # 5. Atomic symlink swap + save watermark state
    atomic_symlink_swap(release_dir, CURRENT_LINK)
    new_state = {
        "tables": {k: (list(v[0]) if isinstance(v[0], tuple) else v[0]) for k, v in compact_queries.items()},
        "samples_by_source": wm["samples_by_source"],
        "artifacts_by_source": wm["artifacts_by_source"],
    }
    save_state(new_state)
    cleanup_old_releases(keep=max(1, keep_releases))
    log.info("Published atomic Parquet release -> %s", CURRENT_LINK)
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description="Atomic, source-partitioned Parquet exporter")
    p.add_argument("--full", action="store_true", help="Force full re-export of all partitions")
    p.add_argument("--dry-run", action="store_true", help="Inspect SQLite watermarks without exporting")
    p.add_argument("--keep-releases", type=int, default=2, help="Number of historical releases to retain")
    args = p.parse_args()
    return run_export(full=args.full, dry_run=args.dry_run, keep_releases=args.keep_releases)


if __name__ == "__main__":
    sys.exit(main())
