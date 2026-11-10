#!/usr/bin/env python3
"""lit_get: everything saved about one paper (record, files, versions, wiki page, optional text)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "lib"))
sys.path.insert(0, str(_HERE))

from common import emit, fail, open_project_db  # noqa: E402
from db import (  # noqa: E402
    derive_review_reasons,
    normalize_arxiv,
    normalize_doi,
    normalize_openalex,
    normalize_pmcid,
    normalize_pmid,
)
from fulltext import local_text, saved_files  # noqa: E402
from wiki import safe_id, vault_root  # noqa: E402


def _candidates(token: str) -> list[str]:
    t = token.strip()
    out = [t]
    for prefix, norm in (
        ("doi", normalize_doi), ("pmid", normalize_pmid), ("pmcid", normalize_pmcid),
        ("arxiv", normalize_arxiv), ("openalex", normalize_openalex),
    ):
        v = norm(t)
        if v:
            out.append(f"{prefix}:{v}")
    return out


def _find(conn, token: str):
    for cand in _candidates(token):
        row = conn.execute("SELECT * FROM publications WHERE id = ?", (cand,)).fetchone()
        if row is None:
            hit = conn.execute("SELECT publication_id FROM aliases WHERE alias = ?", (cand,)).fetchone()
            if hit:
                row = conn.execute("SELECT * FROM publications WHERE id = ?", (hit["publication_id"],)).fetchone()
        if row is not None:
            return dict(row)
    return None


def main() -> int:
    p = argparse.ArgumentParser(description="Show one saved paper")
    p.add_argument("--project_root", required=True)
    p.add_argument("--publication_id", required=True, help="Library id, DOI, PMID, PMCID, arXiv or OpenAlex id")
    p.add_argument("--with_fulltext", action="store_true", help="Include text read from a saved PDF")
    p.add_argument("--fulltext_chars", type=int, default=6000)
    p.add_argument("--with_wiki", action="store_true", help="Include the wiki source page text")
    args = p.parse_args()

    conn, err = open_project_db(args.project_root)
    if err:
        return fail(err)
    try:
        pub = _find(conn, args.publication_id)
        if pub is None:
            return fail(f"no saved paper matches {args.publication_id!r}")
        pid = pub["id"]
        collections = [
            dict(r) for r in conn.execute(
                "SELECT c.id, c.name, c.path FROM collections c "
                "JOIN collection_pubs cp ON cp.collection_id = c.id WHERE cp.publication_id = ?", (pid,),
            )
        ]
        versions = [
            dict(r) for r in conn.execute(
                "SELECT id, from_id, to_id, relation, confidence, needs_confirm FROM publication_versions "
                "WHERE from_id = ? OR to_id = ?", (pid, pid),
            )
        ]
        files = saved_files(conn, pid)
        reasons = derive_review_reasons(conn, pub)
        text = local_text(conn, pid, max_chars=max(500, args.fulltext_chars)) if args.with_fulltext else None
    finally:
        conn.close()

    try:
        authors = json.loads(pub.get("authors_json") or "[]")
    except json.JSONDecodeError:
        authors = []
    try:
        tags = json.loads(pub.get("tags_json") or "[]")
    except json.JSONDecodeError:
        tags = []
    wiki = []
    for c in collections:
        page = vault_root(Path(c["path"])) / "wiki" / "sources" / f"{safe_id(pid)}.md"
        if page.is_file():
            entry = {"collection_id": c["id"], "path": str(page)}
            if args.with_wiki:
                entry["text"] = page.read_text(encoding="utf-8")
            wiki.append(entry)
    emit({
        "ok": True,
        "publication": {
            "id": pid, "title": pub.get("title"), "authors": authors, "year": pub.get("year"),
            "venue": pub.get("venue"), "type": pub.get("type"), "abstract": pub.get("abstract"),
            "doi": pub.get("doi"), "pmid": pub.get("pmid"), "pmcid": pub.get("pmcid"),
            "arxiv": pub.get("arxiv"), "openalex": pub.get("openalex"), "tags": tags,
            "needs_review": bool(pub.get("needs_review")), "review_reasons": reasons,
        },
        "collections": [{"id": c["id"], "name": c["name"]} for c in collections],
        "files": files,
        "versions": versions,
        "wiki_pages": wiki,
        **({"fulltext": text} if text is not None else {}),
    })
    return 0


if __name__ == "__main__":
    sys.exit(main())
