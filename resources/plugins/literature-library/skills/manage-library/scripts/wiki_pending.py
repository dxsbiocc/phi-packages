#!/usr/bin/env python3
"""lit_wiki_pending: list papers that have no wiki page yet and write one compile batch.

The agent writes the pages' prose itself (summary, key claims, concepts) from the batch,
then calls lit_wiki_apply once. One request file covers up to --limit papers.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "lib"))
sys.path.insert(0, str(_HERE))

from common import emit, fail, open_project_db, resolve_collection  # noqa: E402
from fulltext import local_text  # noqa: E402
from wiki import (  # noqa: E402
    INSTRUCTIONS,
    RESPONSE_SCHEMA,
    concepts_needing_definition,
    init_vault,
    request_item,
    safe_id,
    sources_written,
    utc_now,
    vault_root,
)

MAX_BATCH = 25
REQUEST_NAME = "compile_batch.json"
RESPONSE_NAME = "compile_batch_response.json"


def main() -> int:
    p = argparse.ArgumentParser(description="Write a batch of papers for the agent to compile into wiki pages")
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    p.add_argument("--limit", type=int, default=10, help=f"Papers per batch (1-{MAX_BATCH})")
    p.add_argument("--publication_ids", action="append", default=[], help="Rewrite these papers even if a page exists")
    p.add_argument("--with_fulltext", action="store_true", help="Include text from a saved PDF when there is one")
    p.add_argument("--fulltext_chars", type=int, default=6000)
    args = p.parse_args()

    conn, err = open_project_db(args.project_root)
    if err:
        return fail(err)
    try:
        col, err = resolve_collection(conn, args.collection_id or None)
        if err:
            return fail(err)
        collection_id = col["id"]
        col_dir = Path(col["path"])
        init_vault(col_dir, collection_id=collection_id, collection_name=col["name"])
        vault = vault_root(col_dir)

        rows = conn.execute(
            "SELECT p.* FROM publications p JOIN collection_pubs cp ON cp.publication_id = p.id "
            "WHERE cp.collection_id = ? ORDER BY p.year DESC, p.id",
            (collection_id,),
        ).fetchall()
        forced = set(args.publication_ids)
        done = sources_written(vault)

        todo, no_text = [], []
        for r in rows:
            pub = dict(r)
            if pub["id"] in forced or safe_id(pub["id"]) not in done:
                if forced and pub["id"] not in forced:
                    continue
                todo.append(pub)
        limit = max(1, min(int(args.limit), MAX_BATCH))
        items = []
        for pub in todo:
            excerpt = ""
            if args.with_fulltext:
                excerpt = local_text(conn, pub["id"], max_chars=max(500, args.fulltext_chars))["text"]
            if not (pub.get("abstract") or "").strip() and not excerpt:
                no_text.append(pub["id"])
                continue
            if len(items) < limit:
                items.append(request_item(pub, fulltext_excerpt=excerpt))
    finally:
        conn.close()

    request_path = vault / "raw" / REQUEST_NAME
    response_path = vault / "raw" / RESPONSE_NAME
    if items:
        request_path.write_text(
            json.dumps(
                {
                    "schema": "wiki-compile-batch-v1",
                    "created_at": utc_now(),
                    "collection_id": collection_id,
                    "instructions": INSTRUCTIONS,
                    "response_path": str(response_path),
                    "response_schema": RESPONSE_SCHEMA,
                    "items": items,
                },
                ensure_ascii=False,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    n_ready = len(todo) - len(no_text)
    emit(
        {
            "ok": True,
            "collection_id": collection_id,
            "vault_root": str(vault),
            "n_pending": n_ready,
            "n_in_batch": len(items),
            "n_without_text": len(no_text),
            "without_text_ids": no_text[:20],
            "request_path": str(request_path) if items else None,
            "response_path": str(response_path) if items else None,
            "concepts_needing_definition": concepts_needing_definition(vault),
            "next": (
                "Read request_path, write results to response_path, then call lit_wiki_apply."
                if items
                else "Nothing to compile."
                + (" Some papers have no abstract; run lit_curate with action=enrich." if no_text else "")
            ),
        }
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
