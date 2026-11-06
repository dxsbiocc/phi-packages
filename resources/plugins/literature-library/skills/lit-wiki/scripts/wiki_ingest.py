#!/usr/bin/env python3
"""Ingest a publication into the LLM Wiki vault (default compile = agent, v0.4.0).

Agent flow (no plugin API keys):
  1) wiki_ingest ... --compile agent
       → writes compile_request.json, exit 2 (waiting for agent model)
  2) Agent fills compile_response.json with its own model
  3) wiki_ingest ... --compile agent --llm-json compile_response.json
       → validate + write wiki pages
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parent
_DB = Path(__file__).resolve().parents[2] / "literature-library" / "scripts"
sys.path.insert(0, str(_LIB))
sys.path.insert(0, str(_DB))
from wiki_lib import (  # noqa: E402
    collection_dir,
    ingest_publication_to_wiki,
    init_vault,
    safe_id,
    tags_are_empty,
    vault_root,
)
from wiki_llm import (  # noqa: E402
    CompileRequestNeeded,
    build_compile_request,
    resolve_compile_mode,
    write_compile_request,
)
from db import open_db, utc_now, _fts_upsert  # noqa: E402


def _persist_tags(conn, publication_id: str, tags: list[str]) -> bool:
    if not tags:
        return False
    row = conn.execute(
        "SELECT tags_json, title, abstract FROM publications WHERE id = ?",
        (publication_id,),
    ).fetchone()
    if not row:
        return False
    if not tags_are_empty(row["tags_json"]):
        return False
    tags_json = json.dumps(tags, ensure_ascii=False)
    conn.execute(
        "UPDATE publications SET tags_json = ?, updated_at = ? WHERE id = ?",
        (tags_json, utc_now(), publication_id),
    )
    _fts_upsert(conn, publication_id, row["title"] or "", row["abstract"] or "", tags_json)
    conn.commit()
    return True


def _default_request_paths(vault: Path, pub: dict) -> tuple[Path, Path]:
    safe = safe_id(str(pub.get("id") or "unknown"))
    base = vault / "raw" / "publications" / safe
    return base / "compile_request.json", base / "compile_response.json"


def main() -> int:
    p = argparse.ArgumentParser(
        description="Wiki ingest — agent-model compile by default (no plugin API keys)."
    )
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    p.add_argument("--publication_id", default="", help="Canonical publication id in library.sqlite")
    p.add_argument("--raw_path", default="", help="Optional file to copy into raw/")
    p.add_argument("--ensure_vault", action="store_true", help="Create vault if missing")
    p.add_argument(
        "--no_persist_tags",
        action="store_true",
        help="Do not write suggested tags back into SQLite",
    )
    p.add_argument(
        "--compile",
        default="",
        dest="compile_mode",
        help="agent (default) | mock | heuristic | file. Env: LIT_WIKI_COMPILE",
    )
    p.add_argument(
        "--llm-json",
        default="",
        dest="llm_json",
        help="Path to agent-produced compile JSON (apply step)",
    )
    p.add_argument(
        "--apply-compile",
        default="",
        dest="apply_compile",
        help="Alias of --llm-json: apply agent compile_response.json",
    )
    p.add_argument(
        "--emit-request",
        "--request-out",
        default="",
        dest="emit_request",
        help="Write compile_request.json here (aliases: --request-out). Default: vault raw/.../compile_request.json",
    )
    p.add_argument(
        "--request-only",
        action="store_true",
        help="Only emit compile_request.json; never apply (even if --llm-json set)",
    )
    p.add_argument(
        "--allow-heuristic-fallback",
        action="store_true",
        help="On schema/apply failure, fall back to heuristic extractors",
    )
    args = p.parse_args()

    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    if not db_path.exists():
        print(json.dumps({"ok": False, "error": f"missing {db_path}"}))
        return 1
    if not args.publication_id and not args.raw_path:
        print(json.dumps({"ok": False, "error": "provide --publication_id and/or --raw_path"}))
        return 1

    llm_json = args.llm_json or args.apply_compile or ""

    try:
        mode = resolve_compile_mode(args.compile_mode or None, llm_json=llm_json or None)
    except ValueError as e:
        print(json.dumps({"ok": False, "error": str(e)}))
        return 1

    conn = open_db(db_path)
    try:
        row = conn.execute(
            "SELECT id, name, path FROM collections WHERE id = ?", (args.collection_id,)
        ).fetchone()
        if not row:
            print(json.dumps({"ok": False, "error": f"unknown collection {args.collection_id}"}))
            return 1
        col_dir = Path(row["path"]) if row["path"] else collection_dir(root, args.collection_id, db_path)
        vault = vault_root(col_dir)
        if not vault.exists():
            if args.ensure_vault:
                init_vault(col_dir, collection_id=args.collection_id, collection_name=row["name"] or "")
            else:
                print(json.dumps({"ok": False, "error": f"missing vault {vault}; pass --ensure_vault"}))
                return 1

        pub = None
        if args.publication_id:
            prow = conn.execute(
                "SELECT * FROM publications WHERE id = ?", (args.publication_id,)
            ).fetchone()
            if not prow:
                print(json.dumps({"ok": False, "error": f"unknown publication {args.publication_id}"}))
                return 1
            pub = dict(prow)
        else:
            rp = Path(args.raw_path)
            pub = {
                "id": f"file:{rp.stem}",
                "title": rp.stem,
                "abstract": None,
                "tags_json": "[]",
                "needs_review": 0,
            }

        req_default, resp_default = _default_request_paths(vault, pub)

        # Request-only / agent without JSON → emit request
        if args.request_only or (mode == "agent" and not llm_json):
            req_path = Path(args.emit_request) if args.emit_request else req_default
            resp_path = resp_default
            req_path = Path(req_path)
            resp_path = Path(resp_path)
            # When --request-out points outside vault, place response next to the request file.
            if args.emit_request:
                resp_path = Path(args.emit_request).resolve().with_name("compile_response.json")
            req_path.parent.mkdir(parents=True, exist_ok=True)
            resp_path.parent.mkdir(parents=True, exist_ok=True)
            req = build_compile_request(
                pub,
                request_path=str(req_path),
                response_path=str(resp_path),
            )
            write_compile_request(req_path, req)
            out = {
                "ok": True,
                "action": "compile_request",
                "waiting_for_agent": True,
                "publication_id": pub.get("id"),
                "request_path": str(req_path),
                "response_path": str(resp_path),
                "collection_id": args.collection_id,
                "vault_root": str(vault),
                "instructions": req["instructions"],
                "next_step": (
                    f"Synthesize JSON with your model → write {resp_path} → "
                    f"re-run wiki_ingest with --compile agent --llm-json {resp_path}"
                ),
            }
            print(json.dumps(out, ensure_ascii=False))
            return 2

        raw_file = Path(args.raw_path) if args.raw_path else None
        try:
            result = ingest_publication_to_wiki(
                vault,
                pub,
                raw_file=raw_file,
                compile_mode=mode,
                llm_json=llm_json or None,
                allow_heuristic_fallback=bool(args.allow_heuristic_fallback),
            )
        except CompileRequestNeeded as e:
            # Safety net if ingest path forgot llm_json
            req_path = Path(args.emit_request) if args.emit_request else req_default
            write_compile_request(req_path, e.request)
            e.request["request_path"] = str(req_path)
            print(
                json.dumps(
                    {
                        "ok": True,
                        "action": "compile_request",
                        "waiting_for_agent": True,
                        "request_path": str(req_path),
                        "response_path": e.request.get("response_path"),
                        "next_step": e.request.get("next_step"),
                    },
                    ensure_ascii=False,
                )
            )
            return 2
        except Exception as e:  # noqa: BLE001
            print(
                json.dumps(
                    {
                        "ok": False,
                        "error": str(e),
                        "compile_mode": mode,
                        "hint": (
                            "Default is agent compile: emit request (no --llm-json), "
                            "fill compile_response.json with the invoking agent's model, "
                            "then re-run with --llm-json / --apply-compile. "
                            "Offline: --compile mock. Emergency: --compile heuristic. "
                            "No plugin API key is required or supported."
                        ),
                    },
                    ensure_ascii=False,
                )
            )
            return 1

        tags_persisted = False
        if (
            result.get("ok")
            and args.publication_id
            and not args.no_persist_tags
            and result.get("tags_filled")
        ):
            tags_persisted = _persist_tags(
                conn, args.publication_id, list(result.get("suggested_tags") or [])
            )
        result["mode"] = "live"
        result["collection_id"] = args.collection_id
        result["vault_root"] = str(vault)
        result["tags_persisted"] = tags_persisted
        result["compile_mode_resolved"] = mode
        print(json.dumps(result, ensure_ascii=False))
        return 0 if result.get("ok") else 1
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
