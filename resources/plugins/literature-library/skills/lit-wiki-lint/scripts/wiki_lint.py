#!/usr/bin/env python3
"""Lint LLM Wiki vault: orphans, index gaps, broken links, needs_review banners."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parents[2] / "lit-wiki-lib" / "scripts"
_DB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
sys.path.insert(0, str(_LIB))
sys.path.insert(0, str(_DB))
from wiki_lib import collection_dir, lint_vault, vault_root  # noqa: E402
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
            "SELECT path FROM collections WHERE id = ?", (args.collection_id,)
        ).fetchone()
        if not row:
            print(json.dumps({"ok": False, "error": f"unknown collection {args.collection_id}"}))
            return 1
        col_dir = Path(row["path"]) if row["path"] else collection_dir(root, args.collection_id, db_path)
    finally:
        conn.close()

    vault = vault_root(col_dir)
    if not vault.exists():
        print(json.dumps({"ok": False, "error": f"missing vault {vault}"}))
        return 1
    out = lint_vault(vault)
    out["mode"] = "live"
    out["collection_id"] = args.collection_id
    out["vault_root"] = str(vault)
    print(json.dumps(out, ensure_ascii=False))
    # lint reports issues but still ok=True; exit 0 unless vault missing
    return 0


if __name__ == "__main__":
    sys.exit(main())
