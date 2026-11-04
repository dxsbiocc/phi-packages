#!/usr/bin/env python3
"""Wiki agent workflow helper: missing source pages + pending concept merges.

No API keys. Optional --emit_requests batch-writes compile_request.json stubs
for the invoking agent to fill.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

_LIB = Path(__file__).resolve().parents[2] / "lit-wiki-lib" / "scripts"
_DB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
sys.path.insert(0, str(_LIB))
sys.path.insert(0, str(_DB))

from wiki_lib import (  # noqa: E402
    collection_dir,
    list_concept_contributions,
    list_concept_slugs,
    safe_id,
    vault_root,
)
from db import _is_unenriched_shell, open_db  # noqa: E402

HERE = Path(__file__).resolve().parent


def _pending_merges(vault: Path) -> list[dict[str, Any]]:
    raw = vault / "raw"
    concepts = vault / "wiki" / "concepts"
    pending: dict[str, dict[str, Any]] = {}

    if raw.is_dir():
        for path in sorted(raw.glob("concept_merge_request_*.json")):
            slug = path.name[len("concept_merge_request_") : -len(".json")]
            try:
                obj = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                obj = {}
            pending[slug] = {
                "slug": slug,
                "reason": "merge_stub",
                "stub": str(path.relative_to(vault)),
                "n_contributions": len(list_concept_contributions(vault, slug)),
                "concept_name": (obj or {}).get("concept_name") or slug,
            }

    # Also flag ≥2 contribs even without stub
    contrib_root = concepts / "_contrib" if concepts.is_dir() else None
    if contrib_root and contrib_root.is_dir():
        for d in sorted(contrib_root.iterdir()):
            if not d.is_dir():
                continue
            slug = d.name
            n = len(list_concept_contributions(vault, slug))
            if n >= 2 and slug not in pending:
                # Check if Across papers already agent-merged (heuristic: merge label present)
                page = concepts / f"{slug}.md"
                already = False
                if page.is_file():
                    txt = page.read_text(encoding="utf-8")
                    if "merge: `agent-" in txt or "merge: `mock-" in txt:
                        # still list if stub existed; else skip fully merged
                        already = True
                if already:
                    continue
                pending[slug] = {
                    "slug": slug,
                    "reason": "contribs_ge_2",
                    "n_contributions": n,
                    "concept_name": slug,
                }
            elif n >= 2 and slug in pending:
                pending[slug]["n_contributions"] = n

    return sorted(pending.values(), key=lambda x: x["slug"])


def agent_status(
    *,
    project_root: Path,
    collection_id: str,
    limit: int = 50,
) -> dict[str, Any]:
    db_path = project_root / "library.sqlite"
    if not db_path.exists():
        return {"ok": False, "error": f"missing {db_path}"}

    conn = open_db(db_path)
    try:
        row = conn.execute(
            "SELECT path FROM collections WHERE id=?", (collection_id,)
        ).fetchone()
        if not row:
            return {"ok": False, "error": f"unknown collection {collection_id}"}
        col_dir = Path(row["path"]) if row["path"] else collection_dir(
            project_root, collection_id, db_path
        )
        vault = vault_root(col_dir)

        pubs = conn.execute(
            """
            SELECT p.* FROM publications p
            JOIN collection_pubs cp ON cp.publication_id = p.id
            WHERE cp.collection_id = ?
            ORDER BY p.id
            """,
            (collection_id,),
        ).fetchall()
        pubs = [dict(r) for r in pubs]
    finally:
        conn.close()

    sources_dir = vault / "wiki" / "sources" if vault.exists() else None
    existing_sources = set()
    if sources_dir and sources_dir.is_dir():
        existing_sources = {p.stem for p in sources_dir.glob("*.md")}

    missing_wiki: list[dict[str, Any]] = []
    shells: list[str] = []
    for pub in pubs:
        pid = str(pub.get("id") or "")
        if _is_unenriched_shell(pub):
            shells.append(pid)
            continue
        sid = safe_id(pid)
        if sid not in existing_sources:
            missing_wiki.append(
                {
                    "publication_id": pid,
                    "safe_id": sid,
                    "title": (pub.get("title") or "")[:160],
                    "year": pub.get("year"),
                    "has_abstract": bool((pub.get("abstract") or "").strip()),
                }
            )

    merges = _pending_merges(vault) if vault.exists() else []

    suggested: list[str] = []
    for m in missing_wiki[:limit]:
        suggested.append(
            "python skills/lit-wiki-ingest/scripts/wiki_compile_request.py "
            f"--project_root {project_root} --collection_id {collection_id} "
            f"--publication_id {m['publication_id']!r} --ensure_vault"
        )
    for mg in merges[:limit]:
        suggested.append(
            "python skills/lit-wiki-ingest/scripts/concept_merge_request.py "
            f"--project_root {project_root} --collection_id {collection_id} "
            f"--slug {mg['slug']}"
        )
    if missing_wiki:
        suggested.append(
            "# After filling compile_response.json: "
            "python skills/lit-wiki-ingest/scripts/wiki_compile_apply.py "
            f"--project_root {project_root} --collection_id {collection_id} "
            "--publication_id ID --response path/to/compile_response.json"
        )
    if merges:
        suggested.append(
            "# After filling concept_merge_response.json: "
            "python skills/lit-wiki-ingest/scripts/concept_merge_apply.py "
            f"--project_root {project_root} --collection_id {collection_id} "
            "--slug SLUG --llm-json path/to/response.json"
        )

    return {
        "ok": True,
        "action": "wiki_agent_status",
        "project_root": str(project_root),
        "collection_id": collection_id,
        "vault_root": str(vault) if vault.exists() else None,
        "n_publications": len(pubs),
        "n_unenriched_shells": len(shells),
        "n_missing_wiki_sources": len(missing_wiki),
        "missing_wiki_sources": missing_wiki[:limit],
        "n_pending_concept_merges": len(merges),
        "pending_concept_merges": merges[:limit],
        "n_concept_pages": len(list_concept_slugs(vault)) if vault.exists() else 0,
        "suggested_next": suggested[: max(limit * 2, 10)],
        "note": (
            "Default wiki path is agent compile (request → your model JSON → apply). "
            "Use --compile mock only for offline CI — not for live acceptance."
        ),
    }


def emit_compile_requests(
    *,
    project_root: Path,
    collection_id: str,
    missing: list[dict[str, Any]],
    out_dir: Path,
    limit: int,
) -> list[dict[str, Any]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    for item in missing[:limit]:
        pid = item["publication_id"]
        dest = out_dir / f"compile_request_{safe_id(pid)}.json"
        cmd = [
            sys.executable,
            str(HERE / "wiki_compile_request.py"),
            "--project_root",
            str(project_root),
            "--collection_id",
            collection_id,
            "--publication_id",
            pid,
            "--ensure_vault",
            "--emit_request",
            str(dest),
        ]
        # Also accept --request-out via wiki_ingest; try emit_request first
        proc = subprocess.run(cmd, capture_output=True, text=True)
        entry: dict[str, Any] = {
            "publication_id": pid,
            "out": str(dest),
            "returncode": proc.returncode,
        }
        if proc.returncode != 0:
            # fallback: --request-out
            cmd2 = [
                sys.executable,
                str(HERE / "wiki_ingest.py"),
                "--request-only",
                "--compile",
                "agent",
                "--project_root",
                str(project_root),
                "--collection_id",
                collection_id,
                "--publication_id",
                pid,
                "--ensure_vault",
                "--request-out",
                str(dest),
            ]
            proc2 = subprocess.run(cmd2, capture_output=True, text=True)
            entry["returncode"] = proc2.returncode
            entry["stderr"] = (proc2.stderr or proc.stderr or "")[:400]
            try:
                entry["stdout"] = json.loads(proc2.stdout) if proc2.stdout.strip() else None
            except json.JSONDecodeError:
                entry["stdout_raw"] = (proc2.stdout or "")[:400]
        else:
            try:
                entry["stdout"] = json.loads(proc.stdout) if proc.stdout.strip() else None
            except json.JSONDecodeError:
                entry["stdout_raw"] = (proc.stdout or "")[:400]
        entry["written"] = dest.is_file()
        results.append(entry)
    return results


def main() -> int:
    p = argparse.ArgumentParser(
        description="List enriched pubs missing wiki sources + pending concept merges."
    )
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    p.add_argument("--limit", type=int, default=50)
    p.add_argument(
        "--emit_requests",
        action="store_true",
        help="Batch-write compile_request.json for missing wiki sources (agent fills responses)",
    )
    p.add_argument("--out_dir", default="", help="Directory for --emit_requests outputs")
    args = p.parse_args()

    root = Path(args.project_root).resolve()
    out = agent_status(
        project_root=root,
        collection_id=args.collection_id,
        limit=args.limit,
    )
    if not out.get("ok"):
        print(json.dumps(out, ensure_ascii=False))
        return 1

    if args.emit_requests:
        out_dir = Path(args.out_dir).resolve() if args.out_dir else (
            root / "collections" / args.collection_id / "wiki-requests"
        )
        emitted = emit_compile_requests(
            project_root=root,
            collection_id=args.collection_id,
            missing=out.get("missing_wiki_sources") or [],
            out_dir=out_dir,
            limit=args.limit,
        )
        out["emit_requests"] = {
            "out_dir": str(out_dir),
            "n_emitted": sum(1 for e in emitted if e.get("written")),
            "results": emitted,
        }

    print(json.dumps(out, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
