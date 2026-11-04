#!/usr/bin/env python3
"""Apply theme-tree proposals to wiki index (+ optional theme_nodes). Requires --confirm."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parents[2] / "lit-wiki-lib" / "scripts"
_DB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
sys.path.insert(0, str(_LIB))
sys.path.insert(0, str(_DB))
from wiki_lib import apply_theme_tree, collection_dir, vault_root  # noqa: E402
from db import open_db  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    p.add_argument("--proposals_path", required=True, help="JSON from lit-wiki-suggest-tree")
    p.add_argument("--confirm", action="store_true", help="Required to write")
    args = p.parse_args()

    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    prop_path = Path(args.proposals_path)
    if not prop_path.is_file():
        print(json.dumps({"ok": False, "error": f"missing proposals {prop_path}"}))
        return 1
    try:
        payload = json.loads(prop_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        print(json.dumps({"ok": False, "error": f"invalid JSON: {e}"}))
        return 1
    proposals = payload.get("proposals") if isinstance(payload, dict) else payload
    if not isinstance(proposals, list):
        print(json.dumps({"ok": False, "error": "proposals must be a list"}))
        return 1

    if not db_path.exists():
        print(json.dumps({"ok": False, "error": f"missing {db_path}"}))
        return 1

    conn = open_db(db_path)
    try:
        row = conn.execute(
            "SELECT path FROM collections WHERE id = ?", (args.collection_id,)
        ).fetchone()
        if not row:
            print(json.dumps({"ok": False, "error": f"unknown collection {args.collection_id}"}))
            return 1
        col_dir = Path(row["path"]) if row["path"] else collection_dir(root, args.collection_id, db_path)
        vault = vault_root(col_dir)
        if not vault.exists():
            print(json.dumps({"ok": False, "error": f"missing vault {vault}"}))
            return 1

        out = apply_theme_tree(
            vault,
            proposals,
            confirm=bool(args.confirm),
            conn=conn if args.confirm else None,
            collection_id=args.collection_id if args.confirm else None,
        )
        out["mode"] = "live"
        out["collection_id"] = args.collection_id
        out["vault_root"] = str(vault)
        out["proposals_path"] = str(prop_path.resolve())
        print(json.dumps(out, ensure_ascii=False))
        return 0 if out.get("ok") else 1
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
