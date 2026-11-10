#!/usr/bin/env python3
"""lit_query: search what is already saved in the library (metadata and wiki pages).

This never calls an external source; use lit_search for that.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "lib"))
sys.path.insert(0, str(_HERE))

from common import emit, fail, open_project_db, resolve_collection  # noqa: E402
from db import fts_available, query_publications  # noqa: E402
from wiki import query_wiki_markdown, safe_id, vault_root  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description="Search saved papers and wiki pages")
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    p.add_argument("--query", "--q", default="", dest="query", help="Words to find; several words must all match")
    p.add_argument("--tag", default="")
    p.add_argument("--needs_review", default="", help="0 or 1; empty for both")
    p.add_argument("--limit", type=int, default=30)
    p.add_argument("--no_wiki", action="store_true", help="Skip wiki pages")
    args = p.parse_args()

    conn, err = open_project_db(args.project_root)
    if err:
        return fail(err)
    try:
        col, err = resolve_collection(conn, args.collection_id or None)
        if err:
            return fail(err)
        rows = query_publications(
            conn,
            collection_id=col["id"],
            q=args.query or None,
            tag=args.tag or None,
            needs_review=int(args.needs_review) if args.needs_review != "" else None,
            limit=args.limit,
        )
        backend = "fts5" if fts_available(conn) else "like"
    finally:
        conn.close()

    vault = vault_root(Path(col["path"]))
    pages = set(p.stem for p in (vault / "wiki" / "sources").glob("*.md")) if (vault / "wiki" / "sources").is_dir() else set()
    papers = [
        {
            "id": r["id"],
            "title": r.get("title"),
            "year": r.get("year"),
            "venue": r.get("venue"),
            "doi": r.get("doi"),
            "pmid": r.get("pmid"),
            "needs_review": r.get("needs_review"),
            "tags_json": r.get("tags_json"),
            "has_wiki_page": safe_id(r["id"]) in pages,
            "abstract": (r.get("abstract") or "")[:300],
        }
        for r in rows
    ]
    wiki_hits = []
    if args.query and not args.no_wiki and vault.exists():
        wiki_hits = query_wiki_markdown(vault, args.query, limit=args.limit)
    emit({
        "ok": True,
        "collection_id": col["id"],
        "backend": backend,
        "n": len(papers),
        "papers": papers,
        "wiki_hits": wiki_hits,
    })
    return 0


if __name__ == "__main__":
    sys.exit(main())
