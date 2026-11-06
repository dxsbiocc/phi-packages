#!/usr/bin/env python3
"""Initialize LLM Wiki vault beside a literature collection."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parent
_DB = Path(__file__).resolve().parents[2] / "literature-library" / "scripts"
sys.path.insert(0, str(_LIB))
sys.path.insert(0, str(_DB))
from wiki_lib import collection_dir, init_vault  # noqa: E402
from db import open_db  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    args = p.parse_args()
    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    if not db_path.exists():
        print(json.dumps({"ok": False, "error": f"missing {db_path}"}))
        return 1

    conn = open_db(db_path)
    try:
        row = conn.execute(
            "SELECT id, name, path FROM collections WHERE id = ?", (args.collection_id,)
        ).fetchone()
        if not row:
            print(json.dumps({"ok": False, "error": f"unknown collection {args.collection_id}"}))
            return 1
        name = row["name"]
        col_dir = Path(row["path"]) if row["path"] else collection_dir(root, args.collection_id, db_path)
    finally:
        conn.close()

    col_dir.mkdir(parents=True, exist_ok=True)
    out = init_vault(col_dir, collection_id=args.collection_id, collection_name=name or "")
    out["mode"] = "live"
    out["collection_id"] = args.collection_id
    out["collection_name"] = name
    print(json.dumps(out, ensure_ascii=False))
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
