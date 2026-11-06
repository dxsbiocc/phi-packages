#!/usr/bin/env python3
"""Emit grounded concept-merge request(s) for the invoking agent (no API keys)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parent
_DB = Path(__file__).resolve().parents[2] / "literature-library" / "scripts"
sys.path.insert(0, str(_LIB))
sys.path.insert(0, str(_DB))

from wiki_lib import collection_dir, list_concept_contributions, vault_root  # noqa: E402
from wiki_llm import build_concept_merge_request, write_compile_request  # noqa: E402


def _resolve_vault(args: argparse.Namespace) -> Path:
    if args.vault:
        return Path(args.vault).resolve()
    if not args.project_root or not args.collection_id:
        raise SystemExit(json.dumps({"ok": False, "error": "provide --vault or --project_root + --collection_id"}))
    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    col = collection_dir(root, args.collection_id, db_path if db_path.exists() else None)
    return vault_root(col)


def _slugs_with_multi(vault: Path) -> list[str]:
    contrib_root = vault / "wiki" / "concepts" / "_contrib"
    if not contrib_root.is_dir():
        return []
    out = []
    for d in sorted(contrib_root.iterdir()):
        if not d.is_dir():
            continue
        n = len(list(d.glob("*.json")))
        if n >= 2:
            out.append(d.name)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description="Build concept merge request JSON for agent fill.")
    p.add_argument("--vault", default="", help="Path to wiki-vault")
    p.add_argument("--project_root", default="")
    p.add_argument("--collection_id", default="")
    p.add_argument("--slug", default="", help="Concept slug (omit with --all)")
    p.add_argument("--all", action="store_true", help="Emit for every concept with ≥2 contributions")
    p.add_argument("--emit-request", "--request-out", default="", dest="emit_request")
    p.add_argument("--response-path", default="", help="Suggested response path written into request")
    args = p.parse_args()

    try:
        vault = _resolve_vault(args)
    except SystemExit as e:
        print(str(e))
        return 1
    if not vault.exists():
        print(json.dumps({"ok": False, "error": f"missing vault {vault}"}))
        return 1

    if args.all:
        slugs = _slugs_with_multi(vault)
    elif args.slug:
        slugs = [args.slug.strip()]
    else:
        print(json.dumps({"ok": False, "error": "provide --slug or --all"}))
        return 1

    results = []
    for slug in slugs:
        req_path = (
            Path(args.emit_request)
            if args.emit_request and not args.all
            else vault / "raw" / f"concept_merge_request_{slug}.json"
        )
        if args.all and args.emit_request:
            req_path = Path(args.emit_request) / f"concept_merge_request_{slug}.json"
        resp_path = (
            Path(args.response_path)
            if args.response_path and not args.all
            else req_path.with_name(f"concept_merge_response_{slug}.json")
        )
        try:
            req = build_concept_merge_request(
                vault,
                slug,
                request_path=str(req_path),
                response_path=str(resp_path),
            )
        except ValueError as e:
            results.append({"ok": False, "slug": slug, "error": str(e)})
            continue
        write_compile_request(req_path, req)
        results.append(
            {
                "ok": True,
                "slug": slug,
                "request_path": str(req_path),
                "response_path": str(resp_path),
                "n_contributions": req["n_contributions"],
                "requires_needs_review_banner": req["requires_needs_review_banner"],
            }
        )

    ok = all(r.get("ok") for r in results) and bool(results)
    out = {
        "ok": ok,
        "action": "concept_merge_request",
        "waiting_for_agent": ok,
        "vault_root": str(vault),
        "results": results,
        "n": len(results),
    }
    if len(results) == 1 and results[0].get("ok"):
        out.update({k: results[0][k] for k in ("slug", "request_path", "response_path", "n_contributions")})
        out["instructions"] = (
            "Fill response_path with merge JSON using your model, then concept_merge_apply.py --llm-json …"
        )
    print(json.dumps(out, ensure_ascii=False))
    return 2 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
