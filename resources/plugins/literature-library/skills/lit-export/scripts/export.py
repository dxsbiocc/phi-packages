#!/usr/bin/env python3
"""Export a collection to research-evidence search_log.jsonl format."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
sys.path.insert(0, str(_LIB))
from db import open_db, utc_now  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    p.add_argument("--out_path", required=True)
    p.add_argument("--include_abstract", action="store_true", default=True)
    args = p.parse_args()

    db_path = Path(args.project_root).resolve() / "library.sqlite"
    out_path = Path(args.out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    conn = open_db(db_path)
    try:
        rows = conn.execute(
            """
            SELECT p.*, cp.added_at, cp.source_run_id
            FROM publications p
            JOIN collection_pubs cp ON cp.publication_id = p.id
            WHERE cp.collection_id = ?
            ORDER BY cp.added_at
            """,
            (args.collection_id,),
        ).fetchall()
    finally:
        conn.close()

    n = 0
    with out_path.open("w", encoding="utf-8") as f:
        for r in rows:
            url = None
            if r["pmid"]:
                url = f"https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/"
            elif r["doi"]:
                url = f"https://doi.org/{r['doi']}"
            rec = {
                "exported_at": utc_now(),
                "collection_id": args.collection_id,
                "publication_id": r["id"],
                "doi": r["doi"],
                "pmid": r["pmid"],
                "pmcid": r["pmcid"],
                "title": r["title"],
                "year": r["year"],
                "venue": r["venue"],
                "url": url,
                "source_run_id": r["source_run_id"],
                "needs_review": r["needs_review"],
                "query": None,
                "source": "literature-library",
            }
            if args.include_abstract:
                rec["abstract"] = r["abstract"]
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
            n += 1

    print(
        json.dumps(
            {
                "ok": True,
                "mode": "live",
                "collection_id": args.collection_id,
                "n": n,
                "out_path": str(out_path.resolve()),
                "format": "research-evidence.search_log.jsonl",
                "note": "Consumer: research-evidence literature-search / evidence table; this plugin does not write hypothesis cards.",
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
