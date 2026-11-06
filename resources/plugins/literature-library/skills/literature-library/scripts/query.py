#!/usr/bin/env python3
"""Query the literature library (FTS5 or LIKE) with collection/tag/needs_review filters."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parent
sys.path.insert(0, str(_LIB))
from db import fts_available, open_db, query_publications  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    p.add_argument("--q", "--query", default="", dest="q", help="Full-text / LIKE query on title+abstract+tags (alias: --query)")
    p.add_argument("--tag", default="")
    p.add_argument("--needs_review", default="", help="0|1 or empty")
    p.add_argument("--limit", type=int, default=50)
    args = p.parse_args()

    db_path = Path(args.project_root).resolve() / "library.sqlite"
    if not db_path.exists():
        print(json.dumps({"ok": False, "error": f"missing {db_path}"}))
        return 1

    nr = None
    if args.needs_review != "":
        nr = int(args.needs_review)

    conn = open_db(db_path)
    try:
        use_fts = fts_available(conn)
        rows = query_publications(
            conn,
            collection_id=args.collection_id or None,
            q=args.q or None,
            tag=args.tag or None,
            needs_review=nr,
            limit=args.limit,
        )
    finally:
        conn.close()

    # Compact authors for stdout
    pubs = []
    for r in rows:
        pubs.append(
            {
                "id": r["id"],
                "doi": r.get("doi"),
                "pmid": r.get("pmid"),
                "title": r.get("title"),
                "year": r.get("year"),
                "venue": r.get("venue"),
                "needs_review": r.get("needs_review"),
                "tags_json": r.get("tags_json"),
                "abstract": (r.get("abstract") or "")[:400],
            }
        )

    print(
        json.dumps(
            {
                "ok": True,
                "mode": "live",
                "backend": "fts5" if use_fts else "like",
                "n": len(pubs),
                "publications": pubs,
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
