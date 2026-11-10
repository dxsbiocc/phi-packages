#!/usr/bin/env python3
"""Enrich ID-only publication shells with title/abstract via PubMed / Europe PMC.

Default target: pubs whose title is empty/shell and abstract is empty.
Uses resolve_or_create richest-merge; never invents fields.
Mock mode reads tests/fixtures/enrich_mock_zheng.json (or --mock_fixture).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Optional

from db import (  # noqa: E402
    _is_shell_title,
    mark_enriched,
    merge_publication_cluster,
    open_db,
    resolve_or_create,
    set_publication_oa,
)
from oa_verify import verify_record_oa  # noqa: E402
from providers_http import (  # noqa: E402
    TITLE_MATCH_MIN,
    fetch_by_arxiv,
    fetch_by_doi,
    fetch_by_pmid,
    find_by_title,
    title_similarity,
)

DEFAULT_MOCK = (
    Path(__file__).resolve().parents[4] / "tests" / "fixtures" / "enrich_mock_zheng.json"
)


def needs_enrich(row: dict[str, Any]) -> bool:
    """True when a paper lacks what a registry lookup could supply."""
    if row.get("enriched_at"):
        return False  # already looked up; use force to try again
    abstract = (row.get("abstract") or "").strip()
    has_id = bool(row.get("doi") or row.get("pmid") or row.get("arxiv"))
    title = (row.get("title") or "").strip()
    if _is_shell_title(row) or not title:
        return has_id or not abstract
    if has_id:
        # A local PDF carries no author list; registry records always do.
        has_authors = (row.get("authors_json") or "[]") not in ("[]", "null", "")
        return not abstract or not has_authors
    # No identifier: a real, multi-word title and no abstract can still be matched by title.
    return not abstract and len(title.split()) >= 4


def _load_mock(path: Path) -> dict[str, dict[str, Any]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("mock fixture must be object keyed by pmid:… / doi:…")
    return {str(k): v for k, v in raw.items() if isinstance(v, dict)}


def fetch_metadata(
    *,
    pmid: Optional[str],
    doi: Optional[str],
    arxiv: Optional[str] = None,
    title: Optional[str] = None,
    mock: bool = False,
    mock_fixture: Optional[Path] = None,
) -> tuple[Optional[dict[str, Any]], dict[str, Any]]:
    if mock:
        path = Path(mock_fixture) if mock_fixture else DEFAULT_MOCK
        if not path.is_file():
            return None, {"provider": "mock", "message": f"missing fixture {path}"}
        index = _load_mock(path)
        for key in (
            f"pmid:{pmid}" if pmid else None,
            f"doi:{doi}" if doi else None,
            doi,
            pmid,
        ):
            if key and key in index:
                hit = dict(index[key])
                hit.setdefault("source", "mock")
                return hit, {"provider": "mock", "via": key, "fixture": str(path)}
        return None, {"provider": "mock", "message": f"no mock entry for pmid={pmid} doi={doi}"}

    errors: list[str] = []
    if pmid:
        hit, meta = fetch_by_pmid(str(pmid))
        if hit and (hit.get("title") or hit.get("abstract")):
            return hit, meta
        errors.append(f"pmid:{meta.get('message') or 'miss'}")
    if doi:
        hit, meta = fetch_by_doi(str(doi))
        if hit and (hit.get("title") or hit.get("abstract")):
            return hit, meta
        errors.append(f"doi:{meta.get('message') or 'miss'}")
    if arxiv:
        hit, meta = fetch_by_arxiv(str(arxiv))
        if hit and (hit.get("title") or hit.get("abstract")):
            return hit, meta
        errors.append(f"arxiv:{meta.get('message') or 'miss'}")
    if not (pmid or doi or arxiv) and title:
        # No identifier at all (e.g. a PDF without a DOI): accept only a near-identical title.
        hit, meta = find_by_title(title)
        if hit:
            if hit.get("doi"):
                full, fmeta = fetch_by_doi(hit["doi"])  # richer record (abstract) when indexed
                if full and full.get("title") and title_similarity(title, full["title"]) >= TITLE_MATCH_MIN:
                    hit, meta = full, {**fmeta, **{k: meta[k] for k in ("title_similarity",) if k in meta}, "via": "title_then_doi"}
            return hit, meta
        errors.append(f"title:{meta.get('message') or 'miss'}")
    return None, {"provider": "error", "message": "; ".join(errors) or "no strong id"}


def enrich_publications(
    *,
    db_path: Path,
    collection_id: str,
    publication_id: Optional[str] = None,
    dry_run: bool = False,
    mock: bool = False,
    mock_fixture: Optional[Path] = None,
    force: bool = False,
    verify_oa: bool = False,
    oa_mock_map: Optional[dict[str, dict[str, Any]]] = None,
) -> dict[str, Any]:
    conn = open_db(db_path)
    details: list[dict[str, Any]] = []
    oa_checks: list[dict[str, Any]] = []
    n_enriched = n_failed = n_skipped = 0
    try:
        col = conn.execute(
            "SELECT id FROM collections WHERE id=?", (collection_id,)
        ).fetchone()
        if not col:
            return {"ok": False, "error": f"unknown collection {collection_id}"}

        if publication_id:
            rows = conn.execute(
                "SELECT p.* FROM publications p "
                "JOIN collection_pubs cp ON cp.publication_id=p.id "
                "WHERE cp.collection_id=? AND p.id=?",
                (collection_id, publication_id),
            ).fetchall()
            if not rows:
                # allow enrich by id even if not yet linked (still require pub exists)
                rows = conn.execute(
                    "SELECT * FROM publications WHERE id=?", (publication_id,)
                ).fetchall()
        else:
            rows = conn.execute(
                "SELECT p.* FROM publications p "
                "JOIN collection_pubs cp ON cp.publication_id=p.id "
                "WHERE cp.collection_id=?",
                (collection_id,),
            ).fetchall()

        candidates = []
        for r in rows:
            d = dict(r)
            if force or needs_enrich(d):
                candidates.append(d)
            else:
                n_skipped += 1

        for row in candidates:
            pid = row["id"]
            hit, meta = fetch_metadata(
                pmid=row.get("pmid"),
                doi=row.get("doi"),
                arxiv=row.get("arxiv"),
                title=None if _is_shell_title(row) else row.get("title"),
                mock=mock,
                mock_fixture=mock_fixture,
            )
            if not hit:
                n_failed += 1
                if not dry_run and "no near-identical title" in str(meta.get("message")):
                    mark_enriched(conn, pid)  # a definite "not found", not a network error
                details.append(
                    {
                        "publication_id": pid,
                        "action": "failed",
                        "meta": meta,
                    }
                )
                continue
            draft = {
                "source": hit.get("source") or ("mock" if mock else "enrich"),
                "doi": hit.get("doi") or row.get("doi"),
                "pmid": hit.get("pmid") or row.get("pmid"),
                "pmcid": hit.get("pmcid") or row.get("pmcid"),
                "arxiv": hit.get("arxiv") or row.get("arxiv"),
                "title": hit.get("title"),
                "abstract": hit.get("abstract"),
                "authors": hit.get("authors"),
                "year": hit.get("year"),
                "venue": hit.get("venue") or hit.get("journal"),
                "type": hit.get("type") or "journal-article",
                "url": hit.get("url"),
            }
            entry = {
                "publication_id": pid,
                "action": "dry_run" if dry_run else "enriched",
                "via": meta.get("via") or meta.get("provider"),
                "title": (draft.get("title") or "")[:120],
                "abstract_len": len(draft.get("abstract") or ""),
                "meta": meta,
            }
            if dry_run:
                n_enriched += 1
                details.append(entry)
                continue
            result = resolve_or_create(conn, draft)
            entry["merge_action"] = result.get("action")
            entry["canonical_id"] = result.get("publication_id")
            if meta.get("title_similarity") is not None:
                entry["matched_by_title"] = meta["title_similarity"]
            canonical = result.get("publication_id")
            if not dry_run and canonical and canonical != pid and pid.startswith("tmp:"):
                # The temporary record found its identity: fold it into the canonical one.
                merged = merge_publication_cluster(conn, [pid, canonical], preferred_id=canonical, clear_review=True)
                entry["merged_tmp_into"] = merged.get("publication_id") or canonical
            if not dry_run:
                # Found in a registry: do not look this paper up again just because the
                # registry has no abstract for it. Failed lookups stay retriable.
                for mid in {pid, result.get("publication_id")}:
                    if mid:
                        mark_enriched(conn, mid)
            if verify_oa and not dry_run:
                oa_res = verify_record_oa(draft, mock_map=oa_mock_map)
                set_publication_oa(conn, result["publication_id"], oa_res)
                oa_checks.append({"publication_id": result["publication_id"], **oa_res})
                entry["oa"] = {
                    "oa_status": oa_res.get("oa_status"),
                    "url": oa_res.get("url"),
                    "stored_pdf_bytes": False,
                }
            n_enriched += 1
            details.append(entry)

        if not dry_run:
            conn.commit()
        return {
            "ok": True,
            "mode": "live",
            "collection_id": collection_id,
            "n_candidates": len(candidates),
            "n_enriched": n_enriched,
            "n_failed": n_failed,
            "n_skipped": n_skipped,
            "dry_run": dry_run,
            "mock": mock,
            "verify_oa": verify_oa,
            "oa_checks": oa_checks,
            "details": details,
        }
    finally:
        conn.close()


def main() -> int:
    p = argparse.ArgumentParser(description="Enrich ID-only publication shells")
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    p.add_argument("--publication_id", default="")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--mock", action="store_true", help="Use offline mock fixture")
    p.add_argument("--mock_fixture", default="", help="Path to mock JSON fixture")
    p.add_argument(
        "--force",
        action="store_true",
        help="Enrich even if abstract already present",
    )
    p.add_argument(
        "--verify_oa",
        action="store_true",
        help="Soft-verify candidate OA/PDF URLs after enrich; never store PDF bytes",
    )
    p.add_argument(
        "--oa_mock_fixture",
        default="",
        help="JSON map url→mock response for offline OA verify",
    )
    args = p.parse_args()
    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    if not db_path.exists():
        print(json.dumps({"ok": False, "error": f"missing {db_path}"}))
        return 1
    oa_mock_map = None
    if args.oa_mock_fixture:
        oa_mock_map = json.loads(Path(args.oa_mock_fixture).read_text(encoding="utf-8"))
    out = enrich_publications(
        db_path=db_path,
        collection_id=args.collection_id,
        publication_id=args.publication_id or None,
        dry_run=bool(args.dry_run),
        mock=bool(args.mock),
        mock_fixture=Path(args.mock_fixture) if args.mock_fixture else None,
        force=bool(args.force),
        verify_oa=bool(args.verify_oa),
        oa_mock_map=oa_mock_map,
    )
    print(json.dumps(out, ensure_ascii=False))
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
