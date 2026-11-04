#!/usr/bin/env python3
"""Propose theme tree branches from tags/titles → JSON only (no writes)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parents[2] / "lit-wiki-lib" / "scripts"
_DB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
sys.path.insert(0, str(_LIB))
sys.path.insert(0, str(_DB))
from wiki_lib import suggest_theme_tree  # noqa: E402
from db import open_db, query_publications  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    p.add_argument("--max_branches", type=int, default=12)
    p.add_argument("--out_path", default="", help="Optional path to write proposals JSON")
    args = p.parse_args()

    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    if not db_path.exists():
        print(json.dumps({"ok": False, "error": f"missing {db_path}"}))
        return 1

    conn = open_db(db_path)
    try:
        row = conn.execute(
            "SELECT 1 FROM collections WHERE id = ?", (args.collection_id,)
        ).fetchone()
        if not row:
            print(json.dumps({"ok": False, "error": f"unknown collection {args.collection_id}"}))
            return 1
        pubs = query_publications(conn, collection_id=args.collection_id, limit=500)
    finally:
        conn.close()

    out = suggest_theme_tree(pubs, max_branches=args.max_branches)
    out["mode"] = "live"
    out["wrote"] = False
    out["collection_id"] = args.collection_id
    if args.out_path:
        op = Path(args.out_path)
        op.parent.mkdir(parents=True, exist_ok=True)
        op.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        out["out_path"] = str(op.resolve())

    print(json.dumps(out, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
