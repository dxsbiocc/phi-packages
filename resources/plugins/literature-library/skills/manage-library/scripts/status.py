#!/usr/bin/env python3
"""lit_status: one-call overview of a project's library and what needs doing next."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "lib"))
sys.path.insert(0, str(_HERE))

from common import db_path_for, emit, open_project_db  # noqa: E402
from wiki import concepts_needing_definition, safe_id, sources_written, vault_root  # noqa: E402


def _collection_status(conn, col) -> dict:
    cid = col["id"]
    pubs = [dict(r) for r in conn.execute(
        "SELECT p.id, p.title, p.abstract, p.needs_review FROM publications p "
        "JOIN collection_pubs cp ON cp.publication_id = p.id WHERE cp.collection_id = ?", (cid,))]
    ids = [p["id"] for p in pubs]
    files = {}
    if ids:
        for r in conn.execute(
            "SELECT publication_id, file_path FROM pub_files WHERE publication_id IN (%s)" % ",".join("?" * len(ids)), ids
        ):
            files.setdefault(r["publication_id"], []).append(r["file_path"])
    n_missing_files = sum(1 for fs in files.values() for f in fs if not Path(f).is_file())
    no_abstract = [p for p in pubs if not (p["abstract"] or "").strip()]
    pending_links = conn.execute(
        "SELECT COUNT(*) FROM publication_versions WHERE needs_confirm = 1"
    ).fetchone()[0]
    vault = vault_root(Path(col["path"]))
    done = sources_written(vault) if vault.exists() else set()
    n_wiki_pending = sum(1 for p in pubs if safe_id(p["id"]) not in done and (p["abstract"] or "").strip())
    needs_def = concepts_needing_definition(vault) if vault.exists() else []
    runs = [dict(r) for r in conn.execute(
        "SELECT id, query, started_at, status FROM runs WHERE collection_id = ? AND query NOT LIKE 'expand:%' "
        "ORDER BY started_at DESC LIMIT 3", (cid,))]

    n_review = sum(1 for p in pubs if p["needs_review"])
    next_steps = []
    if not pubs:
        next_steps.append("Empty collection: run lit_search or lit_import.")
    if n_review:
        next_steps.append(f"{n_review} paper(s) need an identity decision: lit_curate action=queue.")
    if no_abstract:
        next_steps.append(f"{len(no_abstract)} paper(s) have no abstract: lit_curate action=enrich.")
    if pending_links:
        next_steps.append(f"{pending_links} preprint/journal link(s) await confirmation: lit_curate action=versions.")
    if n_missing_files:
        next_steps.append(f"{n_missing_files} saved file path(s) no longer exist: lit_curate action=rebind.")
    if n_wiki_pending:
        next_steps.append(f"{n_wiki_pending} paper(s) have no wiki page: lit_wiki_pending.")
    if needs_def:
        next_steps.append(f"{len(needs_def)} concept page(s) need a Definition (two or more papers).")
    return {
        "collection_id": cid, "name": col["name"],
        "n_papers": len(pubs), "n_needs_review": n_review, "n_without_abstract": len(no_abstract),
        "n_with_files": len(files), "n_missing_files": n_missing_files,
        "n_pending_version_links": pending_links,
        "wiki": {"n_pages": len(done), "n_pending": n_wiki_pending, "n_concepts_needing_definition": len(needs_def)},
        "recent_searches": runs,
        "next_steps": next_steps,
    }


def main() -> int:
    p = argparse.ArgumentParser(description="Library overview")
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    args = p.parse_args()

    if not db_path_for(args.project_root).exists():
        emit({"ok": True, "initialized": False, "collections": [],
              "next_steps": ["No library yet: create a collection with lit_collection."]})
        return 0
    conn, err = open_project_db(args.project_root)
    try:
        rows = conn.execute("SELECT id, name, path FROM collections ORDER BY created_at").fetchall()
        if args.collection_id:
            rows = [r for r in rows if r["id"] == args.collection_id]
        out = [_collection_status(conn, r) for r in rows]
    finally:
        conn.close()
    emit({"ok": True, "initialized": True, "collections": out,
          "next_steps": [] if out else ["No collection yet: create one with lit_collection."]})
    return 0


if __name__ == "__main__":
    sys.exit(main())
