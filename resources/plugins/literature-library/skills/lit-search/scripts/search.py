#!/usr/bin/env python3
"""Run a literature search for a configured run; write hits.jsonl; auto-ingest.

Default: in-plugin HTTP adapters in providers_http.py (no MCP server in package).
Optional MCP *client*: --backend mcp or LIT_SEARCH_BACKEND=mcp → providers_mcp
calling phi-literature search_publications (stdio bridge / mock fixture).
Offline: --mock. Soft-fallbacks to HTTP when MCP unavailable unless --backend mcp
with force (no auto fallback when user forced mcp and it fails).

0.5.0: freezes executed protocol into run.json; durable hits.jsonl + hits.meta.json;
optional coverage gate flags (or use lit-run-close).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_LIB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
_INGEST = Path(__file__).resolve().parents[2] / "lit-ingest" / "scripts"
sys.path.insert(0, str(_LIB))
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_INGEST))

from db import (  # noqa: E402
    RUN_STATUS_SEARCHED,
    close_run,
    open_db,
    set_run_status,
    utc_now,
)
from providers_http import normalize_source, search_source  # noqa: E402
from providers_mcp import (  # noqa: E402
    mcp_available,
    resolve_backend,
    search_publications_mcp,
)
from ingest import ingest_hits_file  # noqa: E402


def _load_run_config(run_json: Path) -> dict:
    return json.loads(run_json.read_text(encoding="utf-8"))


def _mock_hits(query: str, sources: list[str]) -> list[dict]:
    return [
        {
            "source": sources[0] if sources else "mock",
            "doi": "10.1000/mock.lit.001",
            "pmid": "9990001",
            "title": f"Mock paper matching: {query}",
            "abstract": "Offline mock hit for lit-search smoke tests.",
            "authors": ["Doe, Jane", "Smith, John"],
            "year": 2020,
            "venue": "Mock Journal",
            "type": "journal-article",
            "url": "https://example.org/mock/10.1000/mock.lit.001",
        },
        {
            "source": sources[0] if sources else "mock",
            "doi": "10.1000/mock.lit.002",
            "pmid": "9990002",
            "title": "Another mock exhaustion T cell paper",
            "abstract": "Second mock record.",
            "authors": ["Lee, A"],
            "year": 2021,
            "venue": "Mock Journal",
            "type": "journal-article",
            "url": "https://example.org/mock/10.1000/mock.lit.002",
        },
    ]


def _parse_sources(raw: str | None) -> list[str] | None:
    if not raw:
        return None
    parts = [normalize_source(p) for p in raw.replace(";", ",").split(",") if p.strip()]
    return parts or None


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    p.add_argument("--run_id", required=True)
    p.add_argument("--mock", action="store_true", help="Offline mock hits (no network)")
    p.add_argument(
        "--backend",
        default="",
        help="http (default) | mcp | auto. Or set LIT_SEARCH_BACKEND. Plugin does not ship MCP server.",
    )
    p.add_argument(
        "--mcp_mock_fixture",
        default="",
        help="JSON fixture for offline MCP client tests (search_hits / related_hits)",
    )
    p.add_argument(
        "--mcp_cmd",
        default="",
        help="Override LIT_MCP_CMD stdio JSON-lines bridge command",
    )
    p.add_argument(
        "--mcp_fallback_http",
        action="store_true",
        help="When backend=mcp fails, fall back to HTTP (default for auto; off when mcp forced unless this flag)",
    )
    p.add_argument("--skip_ingest", action="store_true")
    p.add_argument(
        "--sources",
        default=None,
        help="Comma-separated source override. Default: frozen run.json sources.",
    )
    p.add_argument("--min_unique_pubs", type=int, default=None)
    p.add_argument("--max_needs_review", type=int, default=None)
    p.add_argument("--max_version_pending", type=int, default=None)
    p.add_argument("--max_unenriched_shells", type=int, default=None)
    p.add_argument(
        "--fail_under_coverage",
        action="store_true",
        help="Exit 2 when coverage gate breaches (requires a gate threshold)",
    )
    p.add_argument(
        "--close_run",
        action="store_true",
        help="After search+ingest, write coverage.json and mark run status",
    )
    args = p.parse_args()

    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    conn = open_db(db_path)
    try:
        col = conn.execute(
            "SELECT id, path FROM collections WHERE id = ?", (args.collection_id,)
        ).fetchone()
        run = conn.execute("SELECT * FROM runs WHERE id = ?", (args.run_id,)).fetchone()
        if not col or not run:
            print(json.dumps({"ok": False, "error": "collection or run not found"}))
            return 1
        col_path = Path(col["path"])
    finally:
        conn.close()

    run_dir = col_path / "runs" / args.run_id
    run_json = run_dir / "run.json"
    if not run_json.exists():
        print(
            json.dumps(
                {"ok": False, "error": f"missing {run_json}; run lit-run-config first"}
            )
        )
        return 1
    cfg = _load_run_config(run_json)
    query = cfg["query"]
    sources = _parse_sources(args.sources) or [
        normalize_source(s) for s in (cfg.get("sources") or ["pubmed"])
    ]
    limit = int(cfg.get("limit_per_source") or cfg.get("limit") or 20)
    date_from = cfg.get("date_from")
    date_to = cfg.get("date_to")

    all_hits: list[dict] = []
    provider_meta: list[dict] = []
    mode = "live"
    backend = resolve_backend(args.backend or None)
    mcp_banner = (
        "NOTE: Default backend=http (providers_http). "
        "Optional MCP client: --backend mcp / LIT_SEARCH_BACKEND=mcp → "
        "search_publications via LIT_MCP_CMD or --mcp_mock_fixture. "
        "This plugin does not ship an MCP server."
    )

    executed_at = utc_now()
    if args.mock:
        mode = "mock"
        all_hits = _mock_hits(query, sources)
        provider_meta.append({"provider": "mock", "n": len(all_hits), "source": "mock"})
    elif backend == "mcp":
        mock_fx = Path(args.mcp_mock_fixture) if args.mcp_mock_fixture else None
        if not mcp_available(mcp_cmd=args.mcp_cmd or None, mock_fixture=mock_fx):
            if args.mcp_fallback_http or False:
                backend_used = "http_fallback"
                for src in sources:
                    hits, meta = search_source(
                        src, query, limit=limit, date_from=date_from, date_to=date_to
                    )
                    provider_meta.append({**meta, "n": len(hits), "fallback_from": "mcp"})
                    all_hits.extend(hits)
                mode = "partial"
                provider_meta.insert(
                    0,
                    {
                        "provider": "mcp",
                        "ok": False,
                        "error": "mcp_unavailable",
                        "fell_back_to": "http",
                    },
                )
            else:
                print(
                    json.dumps(
                        {
                            "ok": False,
                            "error": "mcp_unavailable",
                            "backend": "mcp",
                            "mcp_banner": mcp_banner,
                            "message": (
                                "MCP backend forced but LIT_MCP_CMD / --mcp_mock_fixture missing. "
                                "Set a bridge, pass --mcp_mock_fixture, use --backend http, "
                                "or --mcp_fallback_http."
                            ),
                        },
                        ensure_ascii=False,
                    )
                )
                return 1
        else:
            hits, meta = search_publications_mcp(
                query,
                sources=sources,
                limit=limit,
                date_from=date_from,
                date_to=date_to,
                mcp_cmd=args.mcp_cmd or None,
                mock_fixture=mock_fx,
            )
            if not meta.get("ok", True) and not hits:
                if args.mcp_fallback_http:
                    for src in sources:
                        h2, m2 = search_source(
                            src, query, limit=limit, date_from=date_from, date_to=date_to
                        )
                        provider_meta.append({**m2, "n": len(h2), "fallback_from": "mcp"})
                        all_hits.extend(h2)
                    mode = "partial"
                    provider_meta.insert(0, {**meta, "fell_back_to": "http"})
                else:
                    print(
                        json.dumps(
                            {
                                "ok": False,
                                "error": meta.get("error") or "mcp_failed",
                                "backend": "mcp",
                                "meta": meta,
                                "mcp_banner": mcp_banner,
                            },
                            ensure_ascii=False,
                        )
                    )
                    return 1
            else:
                mode = "mcp_mock" if meta.get("provider") == "mcp_mock" else "mcp"
                all_hits = hits
                provider_meta.append({**meta, "n": len(hits)})
    else:
        # http (default) or auto→http
        for src in sources:
            hits, meta = search_source(
                src, query, limit=limit, date_from=date_from, date_to=date_to
            )
            provider_meta.append({**meta, "n": len(hits)})
            all_hits.extend(hits)
            if meta.get("provider") == "error" and not hits:
                if mode == "live":
                    mode = "partial"

    # Durable hit log for this run (overwrite = single execution snapshot)
    hits_path = run_dir / "hits.jsonl"
    with hits_path.open("w", encoding="utf-8") as f:
        for h in all_hits:
            f.write(json.dumps(h, ensure_ascii=False) + "\n")
    hits_meta = {
        "run_id": args.run_id,
        "collection_id": args.collection_id,
        "query": query,
        "sources": sources,
        "limit_per_source": limit,
        "date_from": date_from,
        "date_to": date_to,
        "executed_at": executed_at,
        "mode": mode,
        "n_hits": len(all_hits),
        "providers": provider_meta,
        "protocol_version": cfg.get("protocol_version"),
        "frozen_at": cfg.get("frozen_at"),
    }
    (run_dir / "hits.meta.json").write_text(
        json.dumps(hits_meta, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # Update frozen run.json with execution stamp (keep protocol fields)
    cfg_exec = dict(cfg)
    cfg_exec["sources_executed"] = sources
    cfg_exec["executed_at"] = executed_at
    cfg_exec["providers"] = provider_meta
    cfg_exec["n_hits"] = len(all_hits)
    cfg_exec["hits_path"] = str(hits_path)
    cfg_exec["hits_meta_path"] = str(run_dir / "hits.meta.json")
    cfg_exec["status"] = RUN_STATUS_SEARCHED
    cfg_exec["search_mode"] = mode
    run_json.write_text(json.dumps(cfg_exec, ensure_ascii=False, indent=2), encoding="utf-8")

    conn = open_db(db_path)
    try:
        set_run_status(
            conn,
            args.run_id,
            RUN_STATUS_SEARCHED,
            extra_config={
                "sources_executed": sources,
                "executed_at": executed_at,
                "providers": provider_meta,
                "n_hits": len(all_hits),
                "search_mode": mode,
            },
            hits_path=str(hits_path),
        )
        conn.commit()
    finally:
        conn.close()

    ingest_result = None
    if not args.skip_ingest and all_hits:
        ingest_result = ingest_hits_file(
            db_path=db_path,
            collection_id=args.collection_id,
            run_id=args.run_id,
            hits_path=hits_path,
        )

    gate_result = None
    exit_code = 0
    want_close = bool(
        args.close_run
        or args.min_unique_pubs is not None
        or args.max_needs_review is not None
        or args.max_version_pending is not None
        or args.max_unenriched_shells is not None
        or args.fail_under_coverage
    )
    if want_close:
        conn = open_db(db_path)
        try:
            gate_result = close_run(
                conn,
                args.run_id,
                min_unique_pubs=args.min_unique_pubs,
                max_needs_review=args.max_needs_review,
                max_version_pending=args.max_version_pending,
                max_unenriched_shells=args.max_unenriched_shells,
                fail_under_coverage=args.fail_under_coverage,
            )
            conn.commit()
        finally:
            conn.close()
        if gate_result.get("should_fail_cli"):
            exit_code = 2

    out = {
        "ok": True,
        "mode": mode,
        "backend": backend,
        "mcp_banner": mcp_banner,
        "run_id": args.run_id,
        "collection_id": args.collection_id,
        "query": query,
        "sources": sources,
        "n_hits": len(all_hits),
        "hits_path": str(hits_path),
        "hits_meta_path": str(run_dir / "hits.meta.json"),
        "providers": provider_meta,
        "ingest": ingest_result,
        "protocol": {
            "protocol_version": cfg.get("protocol_version"),
            "frozen_at": cfg.get("frozen_at"),
            "executed_at": executed_at,
            "limit_per_source": limit,
            "date_from": date_from,
            "date_to": date_to,
        },
        "gate": (gate_result or {}).get("gate"),
        "run_status": (gate_result or {}).get("status") or RUN_STATUS_SEARCHED,
    }
    print(json.dumps(out, ensure_ascii=False))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
