#!/usr/bin/env python3
"""Repair complementary DOI/PMID (etc.) publication clusters.

Default is dry-run. Pass --apply to write merges.
After apply, wiki source pages are repointed when a vault exists (unless
--no_wiki_repoint).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
# wiki helpers live in the sibling lit-wiki skill
_WIKI = Path(__file__).resolve().parents[2] / "lit-wiki" / "scripts"
sys.path.insert(0, str(_WIKI))
from db import apply_id_merges, drain_id_remaps, find_merge_clusters, open_db, plan_id_merges  # noqa: E402
from wiki_lib import wiki_repoint_remaps  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description="Plan/apply high-confidence strong-ID merges")
    p.add_argument("--project_root", default="", help="Project root containing library.sqlite")
    p.add_argument("--db_path", default="", help="Explicit path to library.sqlite")
    p.add_argument(
        "--apply",
        action="store_true",
        help="Write merges (default: dry-run plan only)",
    )
    p.add_argument(
        "--include_low_confidence",
        action="store_true",
        help="In dry-run, also list non-high-confidence clusters (never applied)",
    )
    p.add_argument(
        "--no_wiki_repoint",
        action="store_true",
        help="Do not rename wiki/sources after id remaps",
    )
    args = p.parse_args()

    if args.db_path:
        db_path = Path(args.db_path)
    elif args.project_root:
        db_path = Path(args.project_root).resolve() / "library.sqlite"
    else:
        print(json.dumps({"ok": False, "error": "provide --project_root or --db_path"}))
        return 1
    if not db_path.exists():
        print(json.dumps({"ok": False, "error": f"missing {db_path}"}))
        return 1

    project_root = Path(args.project_root).resolve() if args.project_root else db_path.parent
    conn = open_db(db_path)
    try:
        drain_id_remaps(conn)
        if args.apply:
            out = apply_id_merges(conn, only_high_confidence=True)
            conn.commit()
            remaps = drain_id_remaps(conn)
            out["mode"] = "apply"
            out["db_path"] = str(db_path.resolve())
            out["id_remaps"] = remaps
            if remaps and not args.no_wiki_repoint:
                out["wiki_repoint"] = wiki_repoint_remaps(
                    project_root, remaps, db_path=db_path, leave_stub=True
                )
            else:
                out["wiki_repoint"] = None
            out["note"] = (
                "Merged complementary strong-ID clusters into DOI-preferred canonical ids. "
                "Wiki sources repointed when vaults existed (stub left at old safe_id)."
            )
        else:
            if args.include_low_confidence:
                clusters = find_merge_clusters(conn)
                out = {
                    "ok": True,
                    "n_clusters": len(clusters),
                    "n_actionable": sum(1 for c in clusters if c.get("high_confidence")),
                    "merges": [
                        {
                            "preferred_id": c["preferred_id"],
                            "absorb": [i for i in c["publication_ids"] if i != c["preferred_id"]],
                            "publication_ids": c["publication_ids"],
                            "high_confidence": c["high_confidence"],
                            "titles": c["titles"],
                        }
                        for c in clusters
                    ],
                }
            else:
                out = plan_id_merges(conn, only_high_confidence=True)
            out["mode"] = "dry-run"
            out["db_path"] = str(db_path.resolve())
            out["note"] = "Dry-run only. Re-run with --apply to write merges."
        print(json.dumps(out, ensure_ascii=False))
        return 0 if out.get("ok") else 1
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
