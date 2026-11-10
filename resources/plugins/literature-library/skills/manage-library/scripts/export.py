#!/usr/bin/env python3
"""lit_export: write the collection as BibTeX, RIS, CSL-JSON, CSV, or a search log."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "lib"))
sys.path.insert(0, str(_HERE))

from common import emit, fail, open_project_db, resolve_collection  # noqa: E402
from db import utc_now  # noqa: E402
from formats import EXTENSIONS, FORMATS, RENDERERS  # noqa: E402

SEARCH_LOG = "search_log"


def _rows(conn, collection_id: str, *, only_ids: list[str], include_review: bool):
    rows = conn.execute(
        "SELECT p.*, cp.added_at, cp.source_run_id FROM publications p "
        "JOIN collection_pubs cp ON cp.publication_id = p.id "
        "WHERE cp.collection_id = ? ORDER BY cp.added_at, p.id",
        (collection_id,),
    ).fetchall()
    out = [dict(r) for r in rows]
    if only_ids:
        out = [r for r in out if r["id"] in set(only_ids)]
    if not include_review:
        out = [r for r in out if not int(r.get("needs_review") or 0)]
    return out


def _search_log(rows: list[dict], collection_id: str) -> str:
    lines = []
    for r in rows:
        url = (
            f"https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/" if r.get("pmid")
            else (f"https://doi.org/{r['doi']}" if r.get("doi") else None)
        )
        lines.append(json.dumps({
            "exported_at": utc_now(), "collection_id": collection_id, "publication_id": r["id"],
            "doi": r.get("doi"), "pmid": r.get("pmid"), "pmcid": r.get("pmcid"), "title": r.get("title"),
            "year": r.get("year"), "venue": r.get("venue"), "url": url,
            "source_run_id": r.get("source_run_id"), "needs_review": r.get("needs_review"),
            "abstract": r.get("abstract"),
        }, ensure_ascii=False))
    return "\n".join(lines) + ("\n" if lines else "")


def main() -> int:
    p = argparse.ArgumentParser(description="Export a collection")
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    p.add_argument("--format", default="bibtex", choices=[*FORMATS, SEARCH_LOG])
    p.add_argument("--out_path", default="", help="Default: <project>/exports/<collection>.<ext>")
    p.add_argument("--publication_ids", action="append", default=[], help="Only these papers")
    p.add_argument("--include_needs_review", action="store_true", help="Include papers with unresolved identity")
    args = p.parse_args()

    conn, err = open_project_db(args.project_root)
    if err:
        return fail(err)
    try:
        col, err = resolve_collection(conn, args.collection_id or None)
        if err:
            return fail(err)
        rows = _rows(conn, col["id"], only_ids=args.publication_ids, include_review=args.include_needs_review)
        n_all = conn.execute(
            "SELECT COUNT(*) FROM collection_pubs WHERE collection_id = ?", (col["id"],)
        ).fetchone()[0]
    finally:
        conn.close()

    ext = ".jsonl" if args.format == SEARCH_LOG else EXTENSIONS[args.format]
    out_path = Path(args.out_path) if args.out_path else (
        Path(args.project_root).resolve() / "exports" / f"{col['id']}{ext}"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    text = _search_log(rows, col["id"]) if args.format == SEARCH_LOG else RENDERERS[args.format](rows)
    out_path.write_text(text, encoding="utf-8")
    emit({
        "ok": True, "collection_id": col["id"], "format": args.format, "n": len(rows),
        "n_skipped_needs_review": 0 if args.include_needs_review or args.publication_ids else n_all - len(rows),
        "out_path": str(out_path.resolve()),
    })
    return 0


if __name__ == "__main__":
    sys.exit(main())
