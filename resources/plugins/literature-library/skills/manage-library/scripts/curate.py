#!/usr/bin/env python3
"""lit_curate: maintenance of what is already saved.

action=queue     list papers that need a human decision (identity conflicts, ID-only shells, pending links)
action=enrich    fill title/abstract of ID-only papers from public metadata
action=versions  find or confirm preprint <-> journal links (see --versions_action)
action=merge_ids merge duplicate records that share complementary strong IDs (dry run unless --apply)
action=rebind    re-attach saved files whose path changed, matched by content hash
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "lib"))
sys.path.insert(0, str(_HERE))

from common import emit, fail, open_project_db, resolve_collection  # noqa: E402

LIB = _HERE / "lib"
ACTIONS = ("queue", "enrich", "versions", "merge_ids", "rebind")


def _run(script: str, argv: list[str]) -> int:
    proc = subprocess.run(
        [sys.executable, str(LIB / script), *argv], capture_output=True, text=True, check=False
    )
    out = proc.stdout.strip()
    try:
        data = json.loads(out.splitlines()[-1]) if out else {}
    except (json.JSONDecodeError, IndexError):
        return fail("tool produced no JSON", stderr=proc.stderr[-500:], stdout=out[-500:])
    if proc.returncode not in (0,) and "ok" not in data:
        data = {"ok": False, "error": proc.stderr[-500:] or "failed"}
    emit(data)
    return 0 if data.get("ok", True) else 1


def main() -> int:
    p = argparse.ArgumentParser(description="Maintain the saved library")
    p.add_argument("--action", required=True, choices=ACTIONS)
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    p.add_argument("--publication_id", default="", help="enrich: only this paper")
    p.add_argument("--limit", type=int, default=100, help="queue: max items")
    p.add_argument("--versions_action", default="detect", choices=("detect", "apply", "propose", "confirm"))
    p.add_argument("--from_id", default="", help="versions propose: preprint id")
    p.add_argument("--to_id", default="", help="versions propose: journal id")
    p.add_argument("--link_id", type=int, default=0, help="versions confirm: link to confirm")
    p.add_argument("--http", action="store_true", help="versions: also ask OpenAlex/Crossref")
    p.add_argument("--auto_merge_high_confidence", action="store_true", help="versions apply: merge sure pairs")
    p.add_argument("--apply", action="store_true", help="merge_ids: write the merges (default: dry run)")
    p.add_argument("--scan_dir", default="", help="rebind: folder holding the moved files")
    p.add_argument("--confirm", action="store_true", help="rebind: write path updates (default: dry run)")
    p.add_argument("--mock", action="store_true", help="Offline fixtures for tests")
    args = p.parse_args()

    conn, err = open_project_db(args.project_root)
    if err:
        return fail(err)
    try:
        col, cerr = resolve_collection(conn, args.collection_id or None)
    finally:
        conn.close()
    cid = col["id"] if col else ""
    if cerr and args.action in ("enrich", "rebind") :
        return fail(cerr)
    root = ["--project_root", args.project_root]
    scope = ["--collection_id", cid] if cid else []

    if args.action == "queue":
        return _run("review_queue.py", [*root, *scope, "--limit", str(args.limit)])
    if args.action == "enrich":
        argv = [*root, *scope]
        if args.publication_id:
            argv += ["--publication_id", args.publication_id]
        if args.mock:
            argv.append("--mock")
        return _run("enrich.py", argv)
    if args.action == "versions":
        argv = [*root, *scope, "--action", args.versions_action]
        for flag, val in (("--from_id", args.from_id), ("--to_id", args.to_id)):
            if val:
                argv += [flag, val]
        if args.link_id:
            argv += ["--link_id", str(args.link_id)]
        for flag, on in (("--http", args.http), ("--auto_merge_high_confidence", args.auto_merge_high_confidence), ("--mock", args.mock)):
            if on:
                argv.append(flag)
        return _run("versions.py", argv)
    if args.action == "merge_ids":
        return _run("merge_ids.py", [*root, *(["--apply"] if args.apply else [])])
    if not args.scan_dir:
        return fail("rebind needs scan_dir")
    argv = [*root, *scope, "--scan_dir", args.scan_dir]
    if args.confirm:
        argv.append("--confirm")
    return _run("rebind_files.py", argv)


if __name__ == "__main__":
    sys.exit(main())
