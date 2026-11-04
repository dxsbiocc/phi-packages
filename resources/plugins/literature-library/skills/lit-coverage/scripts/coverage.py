#!/usr/bin/env python3
"""Coverage report: search hits vs unique ingested pubs + review/shell stats."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_LIB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
sys.path.insert(0, str(_LIB))
from db import coverage_report, open_db  # noqa: E402


def _md_summary(payload: dict) -> str:
    hits = payload.get("hits") or {}
    coll = payload.get("collection") or {}
    revid = payload.get("recall_evidence") or {}
    gold = revid.get("gold_set_recall") or {}
    oow = revid.get("out_of_window") or {}
    unc = revid.get("declared_but_uncontacted_sources") or {}
    lines = [
        "# Literature coverage report",
        "",
        f"- run_id: `{payload.get('run_id') or '—'}`",
        f"- collection_id: `{payload.get('collection_id') or '—'}`",
        f"- evidence_weak: **{payload.get('evidence_weak', gold.get('evidence_weak'))}**",
        "",
        "## Hits",
        f"- n_hits: **{hits.get('n', 0)}**",
        f"- unique_publication_ids: **{hits.get('unique_publication_ids', 0)}**",
        f"- dedupe_action_counts: `{json.dumps(hits.get('dedupe_action_counts') or {}, ensure_ascii=False)}`",
        "",
        "## Collection",
        f"- n_publications: **{coll.get('n_publications', 0)}**",
        f"- n_needs_review: {coll.get('n_needs_review', 0)}",
        f"- n_unenriched_shells: {coll.get('n_unenriched_shells', 0)}",
        f"- n_missing_strong_ids: {coll.get('n_missing_strong_ids', 0)}",
        f"- n_missing_abstract: {coll.get('n_missing_abstract', 0)}",
        f"- n_version_links_pending: {coll.get('n_version_links_pending', 0)}",
        "",
        "## Gold-set recall",
        f"- gold_set_n: **{gold.get('gold_set_n', 0)}**",
        f"- found/missed: **{gold.get('n_found', 0)}** / **{gold.get('n_missed', 0)}**",
        f"- recall: `{gold.get('recall')}`",
        f"- note: {gold.get('note') or '—'}",
        "",
        "## Out-of-window controls",
        f"- n_controls: {oow.get('n_controls', 0)}",
        f"- n_leaked: **{oow.get('n_leaked', 0)}** (filter_ok={oow.get('filter_ok')})",
        f"- note: {oow.get('note') or '—'}",
        "",
        "## Declared-but-uncontacted sources",
        f"- n_uncontacted: **{unc.get('n_uncontacted', 0)}**",
        f"- note: {unc.get('note') or '—'}",
    ]
    for item in unc.get("uncontacted") or []:
        lines.append(f"- `{item.get('source')}`: {item.get('reason')}")
    if not (unc.get("uncontacted") or []):
        lines.append("- (section present; none listed)")
    lines.append("")
    for miss in (gold.get("missed") or [])[:20]:
        lines.append(f"- missed: `{miss.get('id')}`")
    for leak in (oow.get("leaked") or [])[:20]:
        lines.append(f"- leaked OOW: `{leak.get('id')}`")
    if gold.get("missed") or oow.get("leaked"):
        lines.append("")
    hints = payload.get("hints") or []
    if hints:
        lines.append("## Actionable hints")
        for h in hints:
            lines.append(f"- {h}")
        lines.append("")
    return "\n".join(lines)



def main() -> int:
    p = argparse.ArgumentParser(description="Coverage report for a run and/or collection")
    p.add_argument("--project_root", default="")
    p.add_argument("--db_path", default="")
    p.add_argument("--run_id", default="")
    p.add_argument("--collection_id", default="")
    p.add_argument("--md", action="store_true")
    p.add_argument("--md_out", default="")
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
    if not args.run_id and not args.collection_id:
        print(json.dumps({"ok": False, "error": "provide --run_id and/or --collection_id"}))
        return 1

    conn = open_db(db_path)
    try:
        out = coverage_report(
            conn,
            run_id=args.run_id or None,
            collection_id=args.collection_id or None,
        )
        out["db_path"] = str(db_path.resolve())
        out["mode"] = "live"
    finally:
        conn.close()

    if out.get("ok"):
        md = _md_summary(out)
        if args.md_out:
            Path(args.md_out).write_text(md, encoding="utf-8")
            out["md_out"] = str(Path(args.md_out).resolve())
        if args.md:
            print(md, file=sys.stderr)

    print(json.dumps(out, ensure_ascii=False))
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
