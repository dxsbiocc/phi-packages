#!/usr/bin/env python3
"""List publications needing human identity adjudication (not full-text review).

0.5.0: --dig mode is schedulable — empty backlog exits 0 with minimal JSON;
non-empty items include suggested_next_action. Optional --write_report (non-empty only).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
sys.path.insert(0, str(_LIB))
from db import list_human_queue, open_db  # noqa: E402


def _md_summary(payload: dict) -> str:
    lines = [
        f"# Literature review queue",
        "",
        f"- n: **{payload.get('n', 0)}** (total {payload.get('n_total', 0)})",
        f"- collection_id: `{payload.get('collection_id') or '—'}`",
        f"- version_links_pending: {payload.get('version_links_pending', 0)}",
        "",
        "## By reason",
        "",
    ]
    by = payload.get("by_reason") or {}
    if not by:
        lines.append("_empty — backlog=0_")
    else:
        for k, v in sorted(by.items(), key=lambda kv: (-kv[1], kv[0])):
            lines.append(f"- `{k}`: {v}")
    lines.append("")
    items = payload.get("items") or []
    if items:
        lines.append("## Items")
        lines.append("")
        for it in items[:50]:
            reasons = ",".join(it.get("reasons") or [])
            title = (it.get("title") or "")[:80]
            lines.append(
                f"- `{it.get('publication_id')}` · {reasons} · {title}"
            )
            if it.get("suggested_next_action"):
                lines.append(f"  - next: {it['suggested_next_action']}")
            if it.get("related_ids"):
                lines.append(
                    f"  - related: {', '.join('`'+r+'`' for r in it['related_ids'][:8])}"
                )
    return "\n".join(lines) + "\n"


def main() -> int:
    p = argparse.ArgumentParser(description="Human identity review queue")
    p.add_argument("--project_root", default="", help="Project root with library.sqlite")
    p.add_argument("--db_path", default="")
    p.add_argument("--collection_id", default="")
    p.add_argument("--limit", type=int, default=200)
    p.add_argument("--no_shells", action="store_true", help="Skip empty unenriched shells")
    p.add_argument("--no_versions", action="store_true", help="Skip pending version links")
    p.add_argument(
        "--dig",
        action="store_true",
        help="Routine-friendly dig: empty backlog → minimal silent JSON; else items + next actions",
    )
    p.add_argument(
        "--write_report",
        default="",
        help="Write markdown report only when backlog n>0",
    )
    p.add_argument("--md", action="store_true", help="Also print markdown summary to stderr")
    p.add_argument("--md_out", default="", help="Write markdown summary to this path (always)")
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

    conn = open_db(db_path)
    try:
        out = list_human_queue(
            conn,
            collection_id=args.collection_id or None,
            limit=args.limit,
            include_shells=not args.no_shells,
            include_versions=not args.no_versions,
        )
        out["db_path"] = str(db_path.resolve())
        out["mode"] = "dig" if args.dig else "live"
    finally:
        conn.close()

    n = int(out.get("n") or 0)

    if args.dig and n == 0:
        # Silent-friendly for scheduled runs
        silent = {
            "ok": True,
            "dig": True,
            "silent": True,
            "n": 0,
            "n_total": 0,
            "items": [],
            "by_reason": {},
            "backlog_empty": True,
            "collection_id": args.collection_id or None,
            "db_path": out["db_path"],
            "mode": "dig",
        }
        print(json.dumps(silent, ensure_ascii=False))
        return 0

    md = _md_summary(out)
    if args.write_report and n > 0:
        Path(args.write_report).write_text(md, encoding="utf-8")
        out["write_report"] = str(Path(args.write_report).resolve())
    if args.md_out:
        Path(args.md_out).write_text(md, encoding="utf-8")
        out["md_out"] = str(Path(args.md_out).resolve())
    if args.md:
        print(md, file=sys.stderr)

    if args.dig:
        out["dig"] = True
        out["silent"] = False
        out["backlog_empty"] = False

    print(json.dumps(out, ensure_ascii=False))
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
