#!/usr/bin/env python3
"""lit_wiki_apply: validate the agent's batch results and write the wiki pages.

Each result is handled on its own: a malformed one is reported and skipped, the rest are written.
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
from db import _fts_upsert  # noqa: E402
from wiki import (  # noqa: E402
    concepts_needing_definition,
    init_vault,
    parse_tags,
    utc_now,
    validate_compile_payload,
    vault_root,
    write_paper,
)

REQUEST_NAME = "compile_batch.json"
RESPONSE_NAME = "compile_batch_response.json"


def _results(obj) -> list:
    if isinstance(obj, dict) and isinstance(obj.get("results"), list):
        return obj["results"]
    if isinstance(obj, list):
        return obj
    if isinstance(obj, dict):
        return [obj]
    return []


def _find_pub(conn, publication_id: str):
    row = conn.execute("SELECT * FROM publications WHERE id = ?", (publication_id,)).fetchone()
    if row is None:
        alias = conn.execute(
            "SELECT publication_id FROM aliases WHERE alias = ?", (publication_id,)
        ).fetchone()
        if alias:
            row = conn.execute(
                "SELECT * FROM publications WHERE id = ?", (alias["publication_id"],)
            ).fetchone()
    return dict(row) if row else None


def _fill_tags(conn, pub: dict, tags: list[str]) -> bool:
    """Save the agent's tags to SQLite only when the paper has none yet."""
    if not tags or parse_tags(pub.get("tags_json")):
        return False
    tags_json = json.dumps(tags, ensure_ascii=False)
    conn.execute(
        "UPDATE publications SET tags_json = ?, updated_at = ? WHERE id = ?",
        (tags_json, utc_now(), pub["id"]),
    )
    _fts_upsert(conn, pub["id"], pub.get("title") or "", pub.get("abstract") or "", tags_json)
    return True


def main() -> int:
    p = argparse.ArgumentParser(description="Write wiki pages from the agent's compile results")
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    p.add_argument("--response", default="", help="Results JSON; default raw/compile_batch_response.json")
    args = p.parse_args()

    conn, err = open_project_db(args.project_root)
    if err:
        return fail(err)
    try:
        col, err = resolve_collection(conn, args.collection_id or None)
        if err:
            return fail(err)
        col_dir = Path(col["path"])
        init_vault(col_dir, collection_id=col["id"], collection_name=col["name"])
        vault = vault_root(col_dir)

        response = Path(args.response) if args.response else vault / "raw" / RESPONSE_NAME
        if not response.is_file():
            return fail(f"missing {response}; write the results there first")
        try:
            results = _results(json.loads(response.read_text(encoding="utf-8")))
        except json.JSONDecodeError as e:
            return fail(f"response is not valid JSON: {e}")
        if not results:
            return fail("response has no results")

        request_path = vault / "raw" / REQUEST_NAME
        by_id: dict = {}
        if request_path.is_file():
            try:
                for it in json.loads(request_path.read_text(encoding="utf-8")).get("items", []):
                    by_id[it.get("publication_id")] = it
            except json.JSONDecodeError:
                pass

        written, failed = [], []
        for i, raw in enumerate(results):
            pid = str((raw or {}).get("publication_id") or "") if isinstance(raw, dict) else ""
            try:
                payload = validate_compile_payload(raw)
                pub = _find_pub(conn, pid)
                if pub is None:
                    raise ValueError(f"unknown publication_id {pid!r}")
                item = by_id.get(pid) or {}
                basis = item.get("basis", "abstract")
                raw_text = None
                if item.get("fulltext_excerpt"):
                    raw_text = (
                        f"Title: {pub.get('title')}\n\nAbstract:\n{pub.get('abstract') or ''}\n\n"
                        f"Full text excerpt:\n{item['fulltext_excerpt']}\n"
                    )
                res = write_paper(vault, pub, payload, raw_text=raw_text, basis=basis)
                if not res.get("ok"):
                    raise ValueError(res.get("error", "write failed"))
                res["tags_filled"] = _fill_tags(conn, pub, payload["tags"])
                written.append(res)
            except ValueError as e:
                failed.append({"index": i, "publication_id": pid, "error": str(e)})
        conn.commit()
    finally:
        conn.close()

    emit(
        {
            "ok": True,
            "collection_id": col["id"],
            "n_written": len(written),
            "n_failed": len(failed),
            "written": [
                {k: w[k] for k in ("publication_id", "source_page", "concepts_created", "tags_filled")}
                for w in written
            ],
            "failed": failed,
            "concepts_needing_definition": concepts_needing_definition(vault),
            "next": (
                "Fix the failed results and call lit_wiki_apply again."
                if failed
                else "Call lit_wiki_pending for the next batch, or lit_wiki_lint."
            ),
        }
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
