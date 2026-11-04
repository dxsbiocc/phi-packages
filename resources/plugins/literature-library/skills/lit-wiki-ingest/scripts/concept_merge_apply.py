#!/usr/bin/env python3
"""Apply agent/mock concept-merge JSON into wiki/concepts/<slug>.md."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parents[2] / "lit-wiki-lib" / "scripts"
_DB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
sys.path.insert(0, str(_LIB))
sys.path.insert(0, str(_DB))

from wiki_lib import collection_dir, vault_root  # noqa: E402
from wiki_llm import (  # noqa: E402
    apply_concept_merge,
    mock_concept_merge_payload,
)


def _resolve_vault(args: argparse.Namespace) -> Path:
    if args.vault:
        return Path(args.vault).resolve()
    if not args.project_root or not args.collection_id:
        raise SystemExit(json.dumps({"ok": False, "error": "provide --vault or --project_root + --collection_id"}))
    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    col = collection_dir(root, args.collection_id, db_path if db_path.exists() else None)
    return vault_root(col)


def main() -> int:
    p = argparse.ArgumentParser(description="Validate + apply concept merge JSON (or --mock).")
    p.add_argument("--vault", default="")
    p.add_argument("--project_root", default="")
    p.add_argument("--collection_id", default="")
    p.add_argument("--slug", required=True)
    p.add_argument("--llm-json", "--response", default="", dest="llm_json")
    p.add_argument("--mock", action="store_true", help="Deterministic offline merge from contributions")
    args = p.parse_args()

    try:
        vault = _resolve_vault(args)
    except SystemExit as e:
        print(str(e))
        return 1
    if not vault.exists():
        print(json.dumps({"ok": False, "error": f"missing vault {vault}"}))
        return 1

    slug = args.slug.strip()
    try:
        if args.mock:
            payload = mock_concept_merge_payload(vault, slug)
            backend = "mock"
        else:
            if not args.llm_json:
                print(json.dumps({"ok": False, "error": "provide --llm-json or --mock"}))
                return 1
            path = Path(args.llm_json)
            if not path.is_file():
                print(json.dumps({"ok": False, "error": f"llm-json not found: {path}"}))
                return 1
            raw = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(raw, dict) and "definition" not in raw and slug in raw:
                raw = raw[slug]
            payload = raw
            backend = "agent"
        result = apply_concept_merge(vault, slug, payload, backend=backend)
    except Exception as e:  # noqa: BLE001
        print(json.dumps({"ok": False, "error": str(e), "slug": slug}, ensure_ascii=False))
        return 1

    result["vault_root"] = str(vault)
    print(json.dumps(result, ensure_ascii=False))
    return 0 if result.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
