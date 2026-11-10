#!/usr/bin/env python3
"""lit_expand: find papers related to one saved paper (similar / citations / references).

One round per call. By default the hits are only previewed; pass --confirm_write
to save them into the collection (with deduplication). To follow citations further,
call again with a hit from the previous result as the new seed.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import uuid
from pathlib import Path
from typing import Any, Optional

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "lib"))
sys.path.insert(0, str(_HERE))

from common import emit, fail, open_project_db, resolve_collection  # noqa: E402
from db import (  # noqa: E402
    normalize_arxiv,
    normalize_doi,
    normalize_openalex,
    normalize_pmcid,
    normalize_pmid,
    open_db,
    utc_now,
)
from providers_http import dedup_hits, expand_related  # noqa: E402
from ingest import ingest_hits_file  # noqa: E402

MOCK_HITS = [
    {
        "source": "mock",
        "doi": "10.1000/mock.expand.001",
        "pmid": "9991001",
        "title": "Mock related paper (expand similar)",
        "abstract": "Offline mock expand hit.",
        "authors": ["Mock, A"],
        "year": "2018",
        "venue": "Mock Journal",
        "type": "journal-article",
        "url": "https://example.org/mock/expand/001",
    },
    {
        "source": "mock",
        "doi": "10.1000/mock.expand.002",
        "pmid": "9991002",
        "title": "Mock citation / reference expand hit",
        "abstract": "Second offline expand record.",
        "authors": ["Mock, B"],
        "year": "2019",
        "venue": "Mock Journal",
        "type": "journal-article",
        "url": "https://example.org/mock/expand/002",
    },
]


def _clamp_limit(n: int) -> int:
    try:
        v = int(n)
    except (TypeError, ValueError):
        v = 10
    return max(1, min(v, 50))


def _parse_seed_string(seed_id: str) -> dict[str, Optional[str]]:
    """Parse free-form seed into id fields (no network)."""
    s = (seed_id or "").strip()
    out: dict[str, Optional[str]] = {
        "doi": None,
        "pmid": None,
        "pmcid": None,
        "arxiv": None,
        "openalex": None,
        "publication_id": None,
        "raw": s,
    }
    if not s:
        return out

    low = s.lower()
    if low.startswith("doi:"):
        out["doi"] = normalize_doi(s[4:])
        return out
    if low.startswith("pmid:") or low.startswith("pubmed:"):
        out["pmid"] = normalize_pmid(s)
        return out
    if low.startswith("pmcid:") or re.match(r"^pmc\d+$", low):
        out["pmcid"] = normalize_pmcid(s)
        return out
    if low.startswith("arxiv:"):
        out["arxiv"] = normalize_arxiv(s)
        return out
    if low.startswith("openalex:") or re.match(r"^https?://openalex\.org/w\d+", low) or re.match(r"^w\d+$", low):
        out["openalex"] = normalize_openalex(s)
        return out

    doi = normalize_doi(s)
    if doi:
        out["doi"] = doi
        return out
    if re.fullmatch(r"\d{1,9}", s):
        out["pmid"] = s
        return out
    if re.fullmatch(r"PMC\d+", s, re.I):
        out["pmcid"] = normalize_pmcid(s)
        return out
    if re.match(r"^\d{4}\.\d{4,5}(v\d+)?$", s) or re.match(r"^[a-z\-]+/\d+$", s, re.I):
        out["arxiv"] = normalize_arxiv(s)
        return out
    # library publication UUID / canonical id
    out["publication_id"] = s
    return out


def resolve_seed(
    *,
    db_path: Path,
    seed_id: str,
) -> dict[str, Any]:
    """Resolve seed against SQLite when present; else parse string."""
    parsed = _parse_seed_string(seed_id)
    resolved = {
        "seed_id": seed_id,
        "doi": parsed.get("doi"),
        "pmid": parsed.get("pmid"),
        "pmcid": parsed.get("pmcid"),
        "arxiv": parsed.get("arxiv"),
        "openalex": parsed.get("openalex"),
        "publication_id": None,
        "title": None,
        "from_db": False,
    }
    if not db_path.is_file():
        return resolved

    conn = open_db(db_path)
    try:
        row = None
        # Direct publication id
        cand = parsed.get("publication_id") or seed_id
        row = conn.execute("SELECT * FROM publications WHERE id = ?", (cand,)).fetchone()
        if row is None:
            # Try aliases
            alias_candidates: list[str] = []
            if parsed.get("doi"):
                alias_candidates.append(f"doi:{parsed['doi']}")
            if parsed.get("pmid"):
                alias_candidates.append(f"pmid:{parsed['pmid']}")
            if parsed.get("pmcid"):
                alias_candidates.append(f"pmcid:{parsed['pmcid']}")
            if parsed.get("arxiv"):
                alias_candidates.append(f"arxiv:{parsed['arxiv']}")
            if parsed.get("openalex"):
                alias_candidates.append(f"openalex:{parsed['openalex']}")
            alias_candidates.append(seed_id)
            for a in alias_candidates:
                hit = conn.execute(
                    "SELECT publication_id FROM aliases WHERE alias = ?", (a,)
                ).fetchone()
                if hit:
                    row = conn.execute(
                        "SELECT * FROM publications WHERE id = ?", (hit["publication_id"],)
                    ).fetchone()
                    if row:
                        break
        if row:
            d = dict(row)
            resolved["from_db"] = True
            resolved["publication_id"] = d.get("id")
            resolved["doi"] = normalize_doi(d.get("doi")) or resolved["doi"]
            resolved["pmid"] = normalize_pmid(d.get("pmid")) or resolved["pmid"]
            resolved["pmcid"] = normalize_pmcid(d.get("pmcid")) or resolved["pmcid"]
            resolved["arxiv"] = normalize_arxiv(d.get("arxiv")) or resolved["arxiv"]
            resolved["openalex"] = normalize_openalex(d.get("openalex")) or resolved["openalex"]
            resolved["title"] = d.get("title")
    finally:
        conn.close()
    return resolved


def _mock_expand(relation: str, limit: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    hits = dedup_hits(MOCK_HITS, limit=limit)
    for h in hits:
        h = dict(h)
        h["expand_relation"] = relation
    return hits, {
        "ok": True,
        "relation": relation,
        "provider": "mock",
        "via": "fixture",
        "n_candidates": len(hits),
        "providers": [{"provider": "mock", "n": len(hits)}],
    }



def main() -> int:
    p = argparse.ArgumentParser(description="Find papers related to one saved paper")
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    p.add_argument("--seed_id", required=True, help="DOI / PMID / PMCID / arXiv / OpenAlex W... / library id")
    p.add_argument("--relation", default="similar", choices=["similar", "citations", "references"])
    p.add_argument("--limit", type=int, default=10)
    p.add_argument("--confirm_write", action="store_true", help="Save hits into the collection (default: preview)")
    p.add_argument("--mock", action="store_true", help="Offline fake hits for tests")
    p.add_argument("--no_enrich", action="store_true", help="Skip filling metadata of ID-only hits after saving")
    args = p.parse_args()

    limit = _clamp_limit(args.limit)
    db_path = Path(args.project_root).resolve() / "library.sqlite"
    seed = resolve_seed(db_path=db_path, seed_id=args.seed_id)
    has_id = any(seed.get(k) for k in ("doi", "pmid", "openalex", "pmcid", "arxiv"))
    if not has_id and not args.mock:
        return fail("could not resolve seed to DOI/PMID/OpenAlex/PMCID/arXiv", seed=seed)

    if args.mock:
        hits, meta = _mock_expand(args.relation, limit)
    else:
        hits, meta = expand_related(
            relation=args.relation,
            limit=limit,
            pmid=seed.get("pmid"),
            doi=seed.get("doi"),
            openalex=seed.get("openalex"),
            pmcid=seed.get("pmcid"),
            arxiv=seed.get("arxiv"),
        )
    hits = dedup_hits(hits, limit=limit)
    for h in hits:
        h.setdefault("expand_relation", args.relation)
        h.setdefault("expand_seed", args.seed_id)

    out: dict[str, Any] = {
        "ok": bool(meta.get("ok", True)) or bool(hits),
        "seed_id": args.seed_id,
        "seed": seed,
        "relation": args.relation,
        "n_hits": len(hits),
        "hits": hits,
        "wrote": False,
        "providers": meta.get("providers") or [],
    }
    if not hits and meta.get("ok") is False:
        out["ok"] = False
        out["error"] = meta.get("message") or "no related hits"

    if not args.confirm_write:
        out["note"] = "Preview only. Re-run with confirm_write to save these into the collection."
        emit(out)
        return 0

    conn, err = open_project_db(args.project_root)
    if err:
        return fail(err)
    try:
        col, err = resolve_collection(conn, args.collection_id or None)
    finally:
        conn.close()
    if err:
        return fail(err)
    collection_id, col_path = col["id"], Path(col["path"])

    expand_id = f"exp-{utc_now().replace(':', '').replace('-', '')}-{uuid.uuid4().hex[:6]}"
    exp_dir = col_path / "expands" / expand_id
    exp_dir.mkdir(parents=True, exist_ok=True)
    hits_file = exp_dir / "hits.jsonl"
    with hits_file.open("w", encoding="utf-8") as f:
        for h in hits:
            f.write(json.dumps(h, ensure_ascii=False) + "\n")
    (exp_dir / "expand.json").write_text(
        json.dumps(
            {"expand_id": expand_id, "seed": seed, "relation": args.relation, "limit": limit,
             "n_hits": len(hits), "created_at": utc_now()},
            ensure_ascii=False, indent=2,
        ),
        encoding="utf-8",
    )
    out.update(wrote=True, expand_id=expand_id, hits_path=str(hits_file))

    if hits:
        conn = open_db(db_path)
        try:
            # run_hits needs a runs row, so an expansion is stored as a synthetic run.
            conn.execute(
                """
                INSERT OR REPLACE INTO runs(
                  id, collection_id, query, sources_json, config_json,
                  started_at, finished_at, hits_path
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    expand_id, collection_id, f"expand:{args.relation}:{args.seed_id}",
                    json.dumps(["expand", args.relation]),
                    json.dumps({"kind": "expand", "seed_id": args.seed_id,
                                "relation": args.relation, "limit": limit}),
                    utc_now(), utc_now(), str(hits_file),
                ),
            )
            conn.commit()
        finally:
            conn.close()
        out["ingest"] = ingest_hits_file(
            db_path=db_path, collection_id=collection_id, run_id=expand_id, hits_path=hits_file
        )
        if not args.no_enrich:
            try:
                from enrich import enrich_publications  # noqa: E402

                out["enrich"] = enrich_publications(
                    db_path=db_path, collection_id=collection_id, dry_run=False, mock=args.mock
                )
            except Exception as e:  # noqa: BLE001
                out["enrich"] = {"ok": False, "error": str(e)}
    emit(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
