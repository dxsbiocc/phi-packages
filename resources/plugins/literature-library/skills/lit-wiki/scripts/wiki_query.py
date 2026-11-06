#!/usr/bin/env python3
"""Keyword / FTS query over wiki markdown + SQLite metadata. No vectors."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parent
_DB = Path(__file__).resolve().parents[2] / "literature-library" / "scripts"
sys.path.insert(0, str(_LIB))
sys.path.insert(0, str(_DB))
from wiki_lib import collection_dir, query_wiki_markdown, vault_root  # noqa: E402
from db import fts_available, open_db, query_publications  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(
        description="Keyword / FTS query over wiki + SQLite. Multiple terms are AND-matched."
    )
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    p.add_argument(
        "--q",
        "--query",
        required=True,
        dest="q",
        help="Keyword query; space-separated terms use AND semantics (alias: --query)",
    )
    p.add_argument("--limit", type=int, default=50)
    p.add_argument("--wiki_only", action="store_true", help="Skip SQLite metadata query")
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
        vault = vault_root(col_dir)
        wiki_hits = query_wiki_markdown(vault, args.q, limit=args.limit) if vault.exists() else []

        pubs = []
        backend = "none"
        if not args.wiki_only:
            backend = "fts5" if fts_available(conn) else "like"
            pubs = query_publications(
                conn,
                collection_id=args.collection_id,
                q=args.q,
                limit=args.limit,
            )
            pubs = [
                {
                    "id": r["id"],
                    "title": r.get("title"),
                    "year": r.get("year"),
                    "needs_review": r.get("needs_review"),
                    "abstract": (r.get("abstract") or "")[:300],
                }
                for r in pubs
            ]
    finally:
        conn.close()

    print(
        json.dumps(
            {
                "ok": True,
                "mode": "live",
                "embeddings": False,
                "q": args.q,
                "wiki_backend": "keyword",
                "sqlite_backend": backend,
                "n_wiki": len(wiki_hits),
                "n_publications": len(pubs),
                "wiki_hits": wiki_hits,
                "publications": pubs,
                "vault_root": str(vault) if vault.exists() else None,
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
