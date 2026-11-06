#!/usr/bin/env python3
"""Background / polled path-watch that rebinds pub_files by content_hash.

Stdlib only (no watchdog). Polls --scan_dir on --interval; --once for a single
pass. Dry-run by default; --confirm writes. Reuses rebind_pub_files.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Optional

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
_DB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
sys.path.insert(0, str(_DB))

from db import open_db, utc_now  # noqa: E402
from rebind_files import rebind_pub_files  # noqa: E402


def watch_rebind(
    *,
    project_root: Path,
    scan_dir: Path,
    collection_id: Optional[str] = None,
    confirm: bool = False,
    interval: float = 30.0,
    once: bool = False,
    max_passes: Optional[int] = None,
    sample_limit: int = 10,
) -> dict[str, Any]:
    """Poll scan_dir and rebind. Returns last pass status (+ history summary)."""
    db_path = project_root / "library.sqlite"
    if not db_path.exists():
        return {"ok": False, "error": f"missing {db_path}", "action": "watch_rebind"}
    if not scan_dir.is_dir():
        return {
            "ok": False,
            "error": f"scan_dir not a directory: {scan_dir}",
            "action": "watch_rebind",
        }

    passes: list[dict[str, Any]] = []
    n = 0
    started = utc_now()
    while True:
        n += 1
        conn = open_db(db_path)
        try:
            result = rebind_pub_files(
                conn,
                scan_dir=scan_dir,
                collection_id=collection_id,
                confirm=confirm,
                sample_limit=sample_limit,
            )
        finally:
            conn.close()
        result = dict(result)
        result["pass"] = n
        result["watched_at"] = utc_now()
        passes.append(
            {
                "pass": n,
                "n_matched": result.get("n_matched"),
                "n_updated": result.get("n_updated"),
                "n_would_update": result.get("n_would_update"),
                "n_scanned": result.get("n_scanned"),
                "dry_run": result.get("dry_run"),
                "watched_at": result["watched_at"],
            }
        )
        if once or (max_passes is not None and n >= max_passes):
            return {
                "ok": True,
                "action": "watch_rebind",
                "mode": "once" if once else "poll_stopped",
                "dry_run": not confirm,
                "confirm": confirm,
                "interval": interval,
                "n_passes": n,
                "started_at": started,
                "finished_at": utc_now(),
                "project_root": str(project_root.resolve()),
                "db_path": str(db_path),
                "scan_dir": str(scan_dir.resolve()),
                "collection_id": collection_id,
                "last": result,
                "passes": passes,
            }
        time.sleep(max(0.1, float(interval)))


def main() -> int:
    p = argparse.ArgumentParser(
        description="Poll scan_dir and rebind pub_files by content_hash (stdlib; dry-run default)."
    )
    p.add_argument("--project_root", required=True)
    p.add_argument("--scan_dir", required=True)
    p.add_argument("--collection_id", default="")
    p.add_argument("--confirm", action="store_true", help="Write path updates (default: dry-run)")
    p.add_argument(
        "--interval",
        type=float,
        default=30.0,
        help="Seconds between passes when not --once (default 30)",
    )
    p.add_argument("--once", action="store_true", help="Single pass then exit")
    p.add_argument(
        "--max_passes",
        type=int,
        default=None,
        help="Stop after N passes (useful for tests); omit for forever",
    )
    p.add_argument("--sample_limit", type=int, default=10)
    args = p.parse_args()

    root = Path(args.project_root).resolve()
    scan_dir = Path(args.scan_dir).resolve()
    out = watch_rebind(
        project_root=root,
        scan_dir=scan_dir,
        collection_id=args.collection_id or None,
        confirm=bool(args.confirm),
        interval=float(args.interval),
        once=bool(args.once),
        max_passes=args.max_passes,
        sample_limit=args.sample_limit,
    )
    print(json.dumps(out, ensure_ascii=False))
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
