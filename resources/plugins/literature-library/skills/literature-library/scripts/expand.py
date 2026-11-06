#!/usr/bin/env python3
"""Seed → related / citations / references expansion.

Default: HTTP via providers_http. Optional MCP client (--backend mcp) for
phi-literature related_publications. Preview hits + meta JSON; never write
DB until --confirm_write. Plugin does not ship an MCP server.
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
_DB = Path(__file__).resolve().parent
_SEARCH = Path(__file__).resolve().parent
_INGEST = Path(__file__).resolve().parent
sys.path.insert(0, str(_DB))
sys.path.insert(0, str(_SEARCH))
sys.path.insert(0, str(_INGEST))

from db import (  # noqa: E402
    assess_snowball_saturation,
    normalize_arxiv,
    normalize_doi,
    normalize_openalex,
    normalize_pmcid,
    normalize_pmid,
    open_db,
    utc_now,
)
from providers_http import dedup_hits, expand_related  # noqa: E402
from providers_mcp import (  # noqa: E402
    mcp_available,
    related_publications_mcp,
    resolve_backend,
)
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



def _hit_key(h: dict[str, Any]) -> str:
    doi = normalize_doi(h.get("doi"))
    if doi:
        return f"doi:{doi}"
    pmid = normalize_pmid(h.get("pmid"))
    if pmid:
        return f"pmid:{pmid}"
    oa = normalize_openalex(h.get("openalex"))
    if oa:
        return f"openalex:{oa}"
    pmcid = normalize_pmcid(h.get("pmcid"))
    if pmcid:
        return f"pmcid:{pmcid}"
    arxiv = normalize_arxiv(h.get("arxiv"))
    if arxiv:
        return f"arxiv:{arxiv}"
    title = (h.get("title") or "").strip().lower()
    return f"title:{title}" if title else f"row:{id(h)}"


def _seed_token_from_hit(h: dict[str, Any]) -> Optional[str]:
    doi = normalize_doi(h.get("doi"))
    if doi:
        return f"doi:{doi}"
    pmid = normalize_pmid(h.get("pmid"))
    if pmid:
        return f"pmid:{pmid}"
    oa = normalize_openalex(h.get("openalex"))
    if oa:
        return f"openalex:{oa}"
    pmcid = normalize_pmcid(h.get("pmcid"))
    if pmcid:
        return pmcid
    arxiv = normalize_arxiv(h.get("arxiv"))
    if arxiv:
        return f"arxiv:{arxiv}"
    return None


def _expand_once(
    *,
    mock: bool,
    relation: str,
    limit: int,
    seed: dict[str, Any],
    backend: str = "http",
    mcp_cmd: Optional[str] = None,
    mcp_mock_fixture: Optional[Path] = None,
    mcp_fallback_http: bool = False,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if mock:
        return _mock_expand(relation, limit)
    if backend == "mcp":
        if not mcp_available(mcp_cmd=mcp_cmd, mock_fixture=mcp_mock_fixture):
            if mcp_fallback_http:
                hits, meta = expand_related(
                    relation=relation,
                    limit=limit,
                    pmid=seed.get("pmid"),
                    doi=seed.get("doi"),
                    openalex=seed.get("openalex"),
                    pmcid=seed.get("pmcid"),
                    arxiv=seed.get("arxiv"),
                )
                meta = dict(meta)
                meta["fallback_from"] = "mcp"
                return hits, meta
            return [], {
                "ok": False,
                "provider": "mcp",
                "error": "mcp_unavailable",
                "message": "LIT_MCP_CMD / --mcp_mock_fixture missing",
            }
        hits, meta = related_publications_mcp(
            relation=relation,
            limit=limit,
            pmid=seed.get("pmid"),
            doi=seed.get("doi"),
            openalex=seed.get("openalex"),
            pmcid=seed.get("pmcid"),
            arxiv=seed.get("arxiv"),
            mcp_cmd=mcp_cmd,
            mock_fixture=mcp_mock_fixture,
        )
        if (not meta.get("ok", True) and not hits) and mcp_fallback_http:
            hits, meta2 = expand_related(
                relation=relation,
                limit=limit,
                pmid=seed.get("pmid"),
                doi=seed.get("doi"),
                openalex=seed.get("openalex"),
                pmcid=seed.get("pmcid"),
                arxiv=seed.get("arxiv"),
            )
            meta2 = dict(meta2)
            meta2["fallback_from"] = "mcp"
            meta2["mcp_error"] = meta.get("error")
            return hits, meta2
        return hits, meta
    return expand_related(
        relation=relation,
        limit=limit,
        pmid=seed.get("pmid"),
        doi=seed.get("doi"),
        openalex=seed.get("openalex"),
        pmcid=seed.get("pmcid"),
        arxiv=seed.get("arxiv"),
    )



def main() -> int:
    p = argparse.ArgumentParser(description="Expand related literature from a seed (HTTP; no MCP)")
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    p.add_argument("--seed_id", required=True, help="DOI / PMID / PMCID / arxiv / openalex W… / library pub id")
    p.add_argument("--relation", default="similar", choices=["similar", "citations", "references"])
    p.add_argument("--limit", type=int, default=10)
    p.add_argument(
        "--confirm_write",
        action="store_true",
        help="Write hits.jsonl under collection expands/ and ingest (default: preview only)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Force preview-only (same as default without --confirm_write)",
    )
    p.add_argument("--mock", action="store_true", help="Offline mock hits (no network)")
    p.add_argument(
        "--backend",
        default="",
        help="http (default) | mcp. Or LIT_SEARCH_BACKEND. Optional MCP client only.",
    )
    p.add_argument("--mcp_mock_fixture", default="")
    p.add_argument("--mcp_cmd", default="")
    p.add_argument(
        "--mcp_fallback_http",
        action="store_true",
        help="Fall back to HTTP when MCP fails",
    )
    p.add_argument(
        "--max_rounds",
        type=int,
        default=1,
        help="Snowball rounds (default 1). Saturation needs consecutive low-yield or no new seeds.",
    )
    p.add_argument(
        "--low_yield_ratio",
        type=float,
        default=0.05,
        help="New-included / seeds below this ratio counts as low-yield (default 0.05).",
    )
    p.add_argument(
        "--seeds_per_round",
        type=int,
        default=5,
        help="Max seeds expanded per round (default 5).",
    )
    p.add_argument("--skip_ingest", action="store_true", help="With --confirm_write, write jsonl only")
    p.add_argument(
        "--no_enrich",
        action="store_true",
        help="Skip optional enrich of shell hits after ingest (default: enrich if available)",
    )
    args = p.parse_args()

    limit = _clamp_limit(args.limit)
    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    write = bool(args.confirm_write) and not args.dry_run
    backend = resolve_backend(args.backend or None)
    mcp_fx = Path(args.mcp_mock_fixture) if args.mcp_mock_fixture else None

    seed = resolve_seed(db_path=db_path, seed_id=args.seed_id)
    if not any([seed.get("doi"), seed.get("pmid"), seed.get("openalex"), seed.get("pmcid"), seed.get("arxiv")]):
        # still allow mock; live needs at least one resolvable id
        if not args.mock:
            print(
                json.dumps(
                    {
                        "ok": False,
                        "stub": False,
                        "mode": "error",
                        "seed_id": args.seed_id,
                        "relation": args.relation,
        "backend": backend,
                        "limit": limit,
                        "seed": seed,
                        "hits": [],
                        "wrote": False,
                        "error": "could not resolve seed to DOI/PMID/OpenAlex/PMCID/arXiv",
                    },
                    ensure_ascii=False,
                )
            )
            return 1

    max_rounds = max(1, int(args.max_rounds))
    seeds_per_round = max(1, min(int(args.seeds_per_round), 20))
    low_yield_ratio = float(args.low_yield_ratio)

    seen_keys: set[str] = set()
    for k, pref in (
        (seed.get("doi"), "doi"),
        (seed.get("pmid"), "pmid"),
        (seed.get("openalex"), "openalex"),
        (seed.get("pmcid"), "pmcid"),
        (seed.get("arxiv"), "arxiv"),
    ):
        if k:
            seen_keys.add(f"{pref}:{k}")

    seed_queue: list[str] = [args.seed_id]
    all_hits: list[dict[str, Any]] = []
    round_stats: list[dict[str, Any]] = []
    mode = "mock" if args.mock else "live"
    last_meta: dict[str, Any] = {}
    providers_acc: list[dict[str, Any]] = []

    for rnum in range(1, max_rounds + 1):
        round_seeds = seed_queue[:seeds_per_round]
        if not round_seeds:
            round_stats.append(
                {
                    "round": rnum,
                    "n_seeds": 0,
                    "n_candidates": 0,
                    "n_new_included": 0,
                    "seeds": [],
                }
            )
            break
        round_new: list[dict[str, Any]] = []
        n_candidates = 0
        for sid in round_seeds:
            s_resolved = (
                seed if sid == args.seed_id else resolve_seed(db_path=db_path, seed_id=sid)
            )
            hits_r, meta_r = _expand_once(
                mock=args.mock,
                relation=args.relation,
                limit=limit,
                seed=s_resolved,
                backend=backend,
                mcp_cmd=args.mcp_cmd or None,
                mcp_mock_fixture=mcp_fx,
                mcp_fallback_http=bool(args.mcp_fallback_http),
            )
            last_meta = meta_r
            providers_acc.extend(meta_r.get("providers") or [])
            if not args.mock and not hits_r and meta_r.get("ok") is False:
                mode = "partial" if mode == "live" else mode
            for h in hits_r:
                h = dict(h)
                h.setdefault("expand_relation", args.relation)
                h.setdefault("expand_seed", sid)
                h["expand_round"] = rnum
                n_candidates += 1
                key = _hit_key(h)
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                round_new.append(h)
        all_hits.extend(round_new)
        round_stats.append(
            {
                "round": rnum,
                "n_seeds": len(round_seeds),
                "n_candidates": n_candidates,
                "n_new_included": len(round_new),
                "seeds": list(round_seeds),
            }
        )
        next_seeds: list[str] = []
        for h in round_new:
            tok = _seed_token_from_hit(h)
            if tok and tok not in next_seeds:
                next_seeds.append(tok)
        seed_queue = next_seeds
        if len(round_new) == 0:
            break
        sat_partial = assess_snowball_saturation(
            round_stats, max_rounds=max_rounds, low_yield_ratio=low_yield_ratio
        )
        if sat_partial["stop_reason"] in {"consecutive_low_yield", "no_new_eligible_seeds"}:
            break

    saturation = assess_snowball_saturation(
        round_stats, max_rounds=max_rounds, low_yield_ratio=low_yield_ratio
    )
    hits = dedup_hits(all_hits, limit=max(limit * max_rounds, limit))
    for h in hits:
        h.setdefault("expand_relation", args.relation)
        h.setdefault("expand_seed", args.seed_id)
    meta = dict(last_meta or {})
    meta["rounds"] = round_stats
    meta["saturation"] = saturation
    if args.mock:
        mode = "mock"
    elif hits:
        mode = "live" if mode != "partial" else "partial"
    else:
        mode = mode if mode == "partial" else "empty"

    expand_id = f"exp-{utc_now().replace(':', '').replace('-', '')}-{uuid.uuid4().hex[:6]}"
    hits_path: Optional[str] = None
    ingest_result = None
    enrich_result = None
    wrote = False

    # Collection path (needed for write)
    col_path: Optional[Path] = None
    if db_path.is_file():
        conn = open_db(db_path)
        try:
            col = conn.execute(
                "SELECT id, path FROM collections WHERE id = ?", (args.collection_id,)
            ).fetchone()
            if col:
                col_path = Path(col["path"])
        finally:
            conn.close()

    if write:
        if col_path is None:
            print(
                json.dumps(
                    {
                        "ok": False,
                        "stub": False,
                        "mode": mode,
                        "seed_id": args.seed_id,
                        "relation": args.relation,
                        "hits": hits,
                        "wrote": False,
                        "error": "collection not found; cannot confirm_write",
                        "meta": meta,
                    },
                    ensure_ascii=False,
                )
            )
            return 1
        exp_dir = col_path / "expands" / expand_id
        exp_dir.mkdir(parents=True, exist_ok=True)
        path = exp_dir / "hits.jsonl"
        with path.open("w", encoding="utf-8") as f:
            for h in hits:
                f.write(json.dumps(h, ensure_ascii=False) + "\n")
        # sidecar meta
        (exp_dir / "expand.json").write_text(
            json.dumps(
                {
                    "expand_id": expand_id,
                    "seed": seed,
                    "relation": args.relation,
                    "limit": limit,
                    "max_rounds": max_rounds,
                    "n_hits": len(hits),
                    "rounds": round_stats,
                    "saturated": bool(saturation.get("saturated")),
                    "stop_reason": saturation.get("stop_reason"),
                    "saturation": saturation,
                    "meta": meta,
                    "created_at": utc_now(),
                },
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
        hits_path = str(path)
        wrote = True
        if not args.skip_ingest and hits:
            # run_hits FK requires a runs row; treat expand as a synthetic run
            conn = open_db(db_path)
            try:
                conn.execute(
                    """
                    INSERT OR REPLACE INTO runs(
                      id, collection_id, query, sources_json, config_json,
                      started_at, finished_at, hits_path
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        expand_id,
                        args.collection_id,
                        f"expand:{args.relation}:{args.seed_id}",
                        json.dumps(["expand", args.relation]),
                        json.dumps(
                            {
                                "kind": "expand",
                                "seed_id": args.seed_id,
                                "relation": args.relation,
                                "limit": limit,
                            }
                        ),
                        utc_now(),
                        utc_now(),
                        hits_path,
                    ),
                )
                conn.commit()
            finally:
                conn.close()
            ingest_result = ingest_hits_file(
                db_path=db_path,
                collection_id=args.collection_id,
                run_id=expand_id,
                hits_path=path,
            )
            if not args.no_enrich:
                try:
                    _ENRICH = Path(__file__).resolve().parent
                    sys.path.insert(0, str(_ENRICH))
                    from enrich import enrich_publications  # noqa: E402

                    enrich_result = enrich_publications(
                        db_path=db_path,
                        collection_id=args.collection_id,
                        dry_run=False,
                        mock=args.mock,
                    )
                except Exception as e:  # noqa: BLE001
                    enrich_result = {"ok": False, "error": str(e)}

    out = {
        "ok": bool(meta.get("ok", True)) or bool(hits) or args.mock,
        "stub": False,
        "mode": mode if not (write and mode == "live") else "live_write",
        "project_root": str(root),
        "collection_id": args.collection_id,
        "seed_id": args.seed_id,
        "seed": seed,
        "relation": args.relation,
        "limit": limit,
        "max_rounds": max_rounds,
        "hits": hits,
        "n_hits": len(hits),
        "wrote": wrote,
        "confirm_write": write,
        "expand_id": expand_id if wrote else None,
        "hits_path": hits_path,
        "ingest": ingest_result,
        "enrich": enrich_result,
        "providers": providers_acc or meta.get("providers") or [],
        "meta": meta,
        "rounds": round_stats,
        "saturated": bool(saturation.get("saturated")),
        "stop_reason": saturation.get("stop_reason"),
        "saturation": saturation,
        "banner": (
            "Preview only — pass --confirm_write to ingest. "
            + (
                f"Saturated ({saturation.get('stop_reason')})."
                if saturation.get("saturated")
                else f"Closure NOT claimed ({saturation.get('stop_reason')})."
            )
            if not write
            else "Wrote expand hits and attempted ingest. "
            + (
                f"Saturated ({saturation.get('stop_reason')})."
                if saturation.get("saturated")
                else f"Closure NOT claimed ({saturation.get('stop_reason')})."
            )
        ),
    }
    if not hits and not args.mock and meta.get("ok") is False:
        out["ok"] = False
        out["error"] = meta.get("message") or "no related hits"
    # Preview with zero hits is still a successful CLI run (exit 0); hard errors returned above.
    print(json.dumps(out, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
