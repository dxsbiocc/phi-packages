#!/usr/bin/env python3
"""Detect / propose / confirm preprint↔journal version links.

Offline heuristics by default. Optional open HTTP (OpenAlex / Crossref) via
--http; --mock uses fixtures (no network). Live HTTP soft-fails per provider.
0.5.0: --preview_fetch pulls metadata-only shells for missing opposite DOIs
(no PDF); confirm merges to journal DOI preferred.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Optional

_ROOT = Path(__file__).resolve().parents[3]
_LIB = Path(__file__).resolve().parent
_SEARCH = Path(__file__).resolve().parent
_WIKI = Path(__file__).resolve().parents[2] / "lit-wiki" / "scripts"
sys.path.insert(0, str(_LIB))
sys.path.insert(0, str(_SEARCH))
sys.path.insert(0, str(_WIKI))

from db import (  # noqa: E402
    apply_version_detection,
    confirm_version_link,
    detect_version_candidates,
    drain_id_remaps,
    merge_version_proposals,
    open_db,
    preview_fetch_opposite,
    proposals_from_doi_ties,
    propose_version_link,
    unresolved_doi_ties,
)
import providers_http as ph  # noqa: E402
from wiki_lib import wiki_repoint_remaps  # noqa: E402


def _load_mock_fixtures(fixture_dir: Path) -> dict[str, Any]:
    """Load optional mock payloads for OpenAlex/Crossref version discovery."""
    out: dict[str, Any] = {}
    oa = fixture_dir / "openalex_preprint_work.json"
    cr = fixture_dir / "crossref_preprint_work.json"
    pf = fixture_dir / "preview_fetch_dois.json"
    if oa.is_file():
        out["openalex"] = json.loads(oa.read_text(encoding="utf-8"))
    if cr.is_file():
        out["crossref"] = json.loads(cr.read_text(encoding="utf-8"))
    if pf.is_file():
        out["preview_fetch"] = json.loads(pf.read_text(encoding="utf-8"))
    return out


def _fetch_doi_metadata(
    doi: str,
    *,
    mock: bool,
    fixture_dir: Optional[Path],
) -> tuple[Optional[dict[str, Any]], dict[str, Any]]:
    """Metadata-only DOI lookup (no PDF). Mock via preview_fetch_dois.json map."""
    if mock:
        fixtures = _load_mock_fixtures(
            fixture_dir or (_ROOT / "tests" / "fixtures" / "providers")
        )
        mapping = fixtures.get("preview_fetch") or {}
        # mapping: { "10.xxx": {hit...} } or { "dois": { ... } }
        if "dois" in mapping and isinstance(mapping["dois"], dict):
            mapping = mapping["dois"]
        key = doi.lower()
        hit = None
        for k, v in mapping.items():
            if str(k).lower() == key:
                hit = dict(v)
                break
        if hit is None:
            # minimal shell from DOI alone
            hit = {
                "doi": doi,
                "title": f"[preview shell] doi:{doi}",
                "type": "journal-article" if not doi.startswith("10.1101/") else "preprint",
            }
            return hit, {"provider": "mock", "source": "preview_fetch", "via": "minimal"}
        hit.setdefault("doi", doi)
        return hit, {"provider": "mock", "source": "preview_fetch", "via": "fixture"}
    return ph.fetch_by_doi(doi)


def _http_ties_for_db(
    conn,
    *,
    collection_id: Optional[str],
    providers: list[str],
    mock: bool,
    fixture_dir: Optional[Path],
    limit: int = 40,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Collect DOI ties from OpenAlex/Crossref for preprint rows; soft-fail."""
    allowed = None
    if collection_id:
        rows = conn.execute(
            "SELECT publication_id FROM collection_pubs WHERE collection_id=?",
            (collection_id,),
        ).fetchall()
        allowed = {r["publication_id"] for r in rows}

    pubs = [dict(r) for r in conn.execute("SELECT * FROM publications").fetchall()]
    if allowed is not None:
        pubs = [p for p in pubs if p["id"] in allowed]

    # Prefer querying preprint-ish rows (and any with DOI)
    seeds = []
    for p in pubs:
        doi = p.get("doi")
        oa = p.get("openalex")
        if not doi and not oa:
            continue
        seeds.append(p)
        if len(seeds) >= limit:
            break

    all_ties: list[dict[str, Any]] = []
    metas: list[dict[str, Any]] = []

    if mock:
        fixtures = _load_mock_fixtures(
            fixture_dir
            or (_ROOT / "tests" / "fixtures" / "providers")
        )
        # Apply fixture extractors once (offline)
        if "openalex" in providers and "openalex" in fixtures:
            for t in ph.extract_openalex_version_ties(fixtures["openalex"]):
                all_ties.append(t)
            metas.append(
                {
                    "provider": "mock",
                    "source": "openalex",
                    "via": "fixture",
                    "n_ties": len(all_ties),
                }
            )
        if "crossref" in providers and "crossref" in fixtures:
            before = len(all_ties)
            # fixture may be bare message or {message: ...}
            raw = fixtures["crossref"]
            work = raw.get("message") if isinstance(raw, dict) and "message" in raw else raw
            if isinstance(work, dict):
                for t in ph.extract_crossref_version_ties(work):
                    all_ties.append(t)
            metas.append(
                {
                    "provider": "mock",
                    "source": "crossref",
                    "via": "fixture",
                    "n_ties": len(all_ties) - before,
                }
            )
        if not fixtures:
            metas.append(
                {
                    "provider": "mock",
                    "source": "none",
                    "message": "no version fixtures found",
                }
            )
        return all_ties, metas

    # Live: query per seed DOI (soft-fail inside discover_version_ties)
    seen_seed: set[str] = set()
    for p in seeds:
        doi = p.get("doi")
        oa = p.get("openalex")
        seed_key = f"doi:{doi}" if doi else f"openalex:{oa}"
        if seed_key in seen_seed:
            continue
        seen_seed.add(seed_key)
        try:
            ties, m = ph.discover_version_ties(
                doi=doi, openalex_id=oa, providers=providers
            )
            all_ties.extend(ties)
            metas.extend(m)
        except Exception as e:  # noqa: BLE001
            metas.append(
                {
                    "provider": "error",
                    "source": "discover",
                    "message": str(e),
                    "seed": seed_key,
                }
            )
    return all_ties, metas


def _maybe_wiki_repoint(
    *,
    project_root: Optional[Path],
    db_path: Path,
    remaps: list[dict[str, str]],
    collection_id: str,
    enabled: bool,
) -> Optional[dict[str, Any]]:
    if not enabled or not remaps:
        return None
    root = project_root or db_path.parent
    return wiki_repoint_remaps(
        root,
        remaps,
        collection_id=collection_id or None,
        db_path=db_path,
        leave_stub=True,
    )


def main() -> int:
    p = argparse.ArgumentParser(description="Preprint↔journal version chain tools")
    p.add_argument("--project_root", default="")
    p.add_argument("--db_path", default="")
    p.add_argument(
        "--action",
        choices=("detect", "propose", "confirm", "apply"),
        default="detect",
        help="detect=dry list; apply=write links (+ optional auto-merge); "
        "propose=manual link; confirm=confirm link_id",
    )
    p.add_argument("--collection_id", default="")
    p.add_argument("--from_id", default="")
    p.add_argument("--to_id", default="")
    p.add_argument("--relation", default="preprint_of")
    p.add_argument("--confidence", type=float, default=0.5)
    p.add_argument("--link_id", type=int, default=0)
    p.add_argument("--no_merge", action="store_true", help="On confirm, do not merge rows")
    p.add_argument(
        "--auto_merge_high_confidence",
        action="store_true",
        help="With --action apply: merge high-confidence pairs into journal DOI",
    )
    p.add_argument("--evidence_json", default="", help="Optional JSON object string for propose")
    p.add_argument(
        "--http",
        action="store_true",
        help="Also query OpenAlex/Crossref for version ties (soft-fail)",
    )
    p.add_argument(
        "--providers",
        default="openalex,crossref",
        help="Comma list used with --http/--mock (default openalex,crossref)",
    )
    p.add_argument(
        "--mock",
        action="store_true",
        help="Offline fixtures instead of live HTTP (implies version HTTP path)",
    )
    p.add_argument(
        "--fixture_dir",
        default="",
        help="Directory with openalex_preprint_work.json / crossref_preprint_work.json",
    )
    p.add_argument(
        "--no_wiki_repoint",
        action="store_true",
        help="Do not rename wiki sources after id merges",
    )
    p.add_argument(
        "--preview_fetch",
        action="store_true",
        help="When opposite DOI missing from DB: fetch metadata-only shell + pending version link (no PDF)",
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
    providers = [x.strip().lower() for x in args.providers.split(",") if x.strip()]
    use_http = bool(args.http or args.mock)
    fixture_dir = Path(args.fixture_dir) if args.fixture_dir else None

    conn = open_db(db_path)
    drain_id_remaps(conn)  # clear any stale
    try:
        http_proposals: list[dict[str, Any]] = []
        http_meta: list[dict[str, Any]] = []
        preview_results: list[dict[str, Any]] = []
        unresolved: list[dict[str, Any]] = []
        if use_http and args.action in ("detect", "apply"):
            ties, http_meta = _http_ties_for_db(
                conn,
                collection_id=args.collection_id or None,
                providers=providers,
                mock=bool(args.mock),
                fixture_dir=fixture_dir,
            )
            http_proposals = proposals_from_doi_ties(conn, ties)
            unresolved = unresolved_doi_ties(conn, ties)
            if args.preview_fetch and unresolved:
                for u in unresolved:
                    hit, meta = _fetch_doi_metadata(
                        u["missing_doi"],
                        mock=bool(args.mock),
                        fixture_dir=fixture_dir,
                    )
                    if args.action == "detect":
                        preview_results.append(
                            {
                                "ok": True,
                                "dry_run": True,
                                "missing_doi": u["missing_doi"],
                                "present_id": u["present_id"],
                                "would_fetch": hit,
                                "fetch_meta": meta,
                            }
                        )
                    else:
                        # apply path: actually create shell + pending link
                        pr = preview_fetch_opposite(
                            conn,
                            missing_doi=u["missing_doi"],
                            present_id=u["present_id"],
                            relation=u.get("relation") or "preprint_of",
                            confidence=float(u.get("confidence") or 0.5),
                            evidence={
                                **(u.get("evidence") or {}),
                                "http_provider": u.get("provider"),
                                "fetch_meta": meta,
                            },
                            collection_id=args.collection_id or None,
                            hit=hit,
                            missing_role=u.get("missing_role") or "to",
                        )
                        preview_results.append(pr)
                if args.action == "apply" and preview_results:
                    conn.commit()
                    # Refresh http proposals now that opposite DOIs may exist
                    http_proposals = proposals_from_doi_ties(conn, ties)

        if args.action == "detect":
            offline = detect_version_candidates(
                conn, collection_id=args.collection_id or None
            )
            proposals = merge_version_proposals(offline, http_proposals)
            out: dict[str, Any] = {
                "ok": True,
                "mode": "detect",
                "n": len(proposals),
                "n_offline": len(offline),
                "n_http": len(http_proposals),
                "http_enabled": use_http,
                "http_mock": bool(args.mock),
                "http_meta": http_meta,
                "proposals": proposals,
                "unresolved_opposite": unresolved,
                "n_unresolved": len(unresolved),
                "preview_fetch": bool(args.preview_fetch),
                "preview_results": preview_results,
                "note": "Dry list only. Use --action apply to write publication_versions rows. "
                "With --preview_fetch, apply creates metadata-only shells for missing opposite DOIs.",
            }
        elif args.action == "apply":
            out = apply_version_detection(
                conn,
                collection_id=args.collection_id or None,
                auto_merge_high_confidence=args.auto_merge_high_confidence,
                extra_proposals=http_proposals,
            )
            conn.commit()
            remaps = drain_id_remaps(conn)
            out["mode"] = "apply"
            out["http_enabled"] = use_http
            out["http_mock"] = bool(args.mock)
            out["n_http"] = len(http_proposals)
            out["http_meta"] = http_meta
            out["unresolved_opposite"] = unresolved
            out["n_unresolved"] = len(unresolved)
            out["preview_fetch"] = bool(args.preview_fetch)
            out["preview_results"] = preview_results
            out["id_remaps"] = remaps
            out["wiki_repoint"] = _maybe_wiki_repoint(
                project_root=project_root,
                db_path=db_path,
                remaps=remaps,
                collection_id=args.collection_id,
                enabled=not args.no_wiki_repoint,
            )
        elif args.action == "propose":
            if not args.from_id or not args.to_id:
                out = {"ok": False, "error": "propose requires --from_id and --to_id"}
            else:
                evidence = None
                if args.evidence_json:
                    evidence = json.loads(args.evidence_json)
                out = propose_version_link(
                    conn,
                    args.from_id,
                    args.to_id,
                    relation=args.relation,
                    confidence=args.confidence,
                    needs_confirm=1,
                    evidence=evidence,
                )
                conn.commit()
                out["mode"] = "propose"
        elif args.action == "confirm":
            if not args.link_id:
                out = {"ok": False, "error": "confirm requires --link_id"}
            else:
                out = confirm_version_link(
                    conn, args.link_id, merge=not args.no_merge
                )
                conn.commit()
                remaps = drain_id_remaps(conn)
                out["mode"] = "confirm"
                out["id_remaps"] = remaps
                out["wiki_repoint"] = _maybe_wiki_repoint(
                    project_root=project_root,
                    db_path=db_path,
                    remaps=remaps,
                    collection_id=args.collection_id,
                    enabled=not args.no_wiki_repoint,
                )
        else:
            out = {"ok": False, "error": f"unknown action {args.action}"}
        out["db_path"] = str(db_path.resolve())
    finally:
        conn.close()

    print(json.dumps(out, ensure_ascii=False))
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
