#!/usr/bin/env python3
"""Close a search run: write coverage.json(+md), evaluate gate, mark status."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parent
sys.path.insert(0, str(_LIB))
from db import close_run, open_db  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description="Coverage gate + mark run passed/needs_attention")
    p.add_argument("--project_root", default="")
    p.add_argument("--db_path", default="")
    p.add_argument("--run_id", required=True)
    p.add_argument("--min_unique_pubs", type=int, default=None)
    p.add_argument("--max_needs_review", type=int, default=None)
    p.add_argument("--max_version_pending", type=int, default=None)
    p.add_argument("--max_unenriched_shells", type=int, default=None)
    p.add_argument(
        "--fail_under_coverage",
        action="store_true",
        help="Exit 2 when gate breaches",
    )
    p.add_argument("--coverage_json", default="", help="Override coverage.json path")
    p.add_argument("--coverage_md", default="", help="Override coverage.md path")
    args = p.parse_args()

    if args.db_path:
        db_path = Path(args.db_path)
    elif args.project_root:
        db_path = Path(args.project_root).resolve() / "library.sqlite"
    else:
        print(json.dumps({"ok": False, "error": "provide --project_root or --db_path"}))
        return 1
    if not db_path.exists():
        print(json.dumps({"ok": False, "error": f"missing {db_path}"}))
        return 1

    conn = open_db(db_path)
    try:
        out = close_run(
            conn,
            args.run_id,
            min_unique_pubs=args.min_unique_pubs,
            max_needs_review=args.max_needs_review,
            max_version_pending=args.max_version_pending,
            max_unenriched_shells=args.max_unenriched_shells,
            fail_under_coverage=args.fail_under_coverage,
            write_coverage_json=Path(args.coverage_json) if args.coverage_json else None,
            write_coverage_md=Path(args.coverage_md) if args.coverage_md else None,
        )
        out["db_path"] = str(db_path.resolve())
        out["mode"] = "live"
        conn.commit()
    finally:
        conn.close()

    print(json.dumps(out, ensure_ascii=False))
    if not out.get("ok"):
        return 1
    if out.get("should_fail_cli"):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
