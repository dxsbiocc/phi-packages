#!/usr/bin/env python3
"""lit_wiki_lint: check the wiki vault for gaps, broken links, orphans and stale banners."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "lib"))
sys.path.insert(0, str(_HERE))

from common import emit, fail, open_project_db, resolve_collection  # noqa: E402
from wiki import lint_vault, vault_root  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description="Check the wiki vault")
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    args = p.parse_args()

    conn, err = open_project_db(args.project_root)
    if err:
        return fail(err)
    try:
        col, err = resolve_collection(conn, args.collection_id or None)
    finally:
        conn.close()
    if err:
        return fail(err)

    vault = vault_root(Path(col["path"]))
    if not vault.exists():
        return fail(f"missing vault {vault}; call lit_wiki_pending to create it")
    out = lint_vault(vault)
    out.update(collection_id=col["id"], vault_root=str(vault))
    emit(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
