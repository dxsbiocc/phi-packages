#!/usr/bin/env python3
"""lit_search: query external literature sources and save the hits into the library.

One call does what used to take three: it records the search (run), fetches hits
from the chosen sources, writes hits.jsonl, and ingests them with deduplication.
Only external sources are queried here; use lit_query to search what is already saved.
"""
from __future__ import annotations

import argparse
import json
import sys
import uuid
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "lib"))
sys.path.insert(0, str(_HERE))

from common import emit, fail, open_project_db, resolve_collection  # noqa: E402
from db import (  # noqa: E402
    RUN_STATUS_SEARCHED,
    build_frozen_protocol,
    open_db,
    set_run_status,
    utc_now,
)
from ingest import ingest_hits_file  # noqa: E402
from providers_http import SUPPORTED_SOURCES, normalize_source, search_source  # noqa: E402

DEFAULT_SOURCES = ("pubmed", "europepmc")
MAX_SOURCES = 4
MAX_LIMIT = 50


def _mock_hits(query: str, source: str) -> list[dict]:
    return [
        {
            "source": source,
            "doi": f"10.1000/mock.lit.00{i}",
            "pmid": f"999000{i}",
            "title": f"Mock paper {i} matching: {query}",
            "abstract": "Offline mock hit for tests.",
            "authors": ["Doe, Jane"],
            "year": 2020 + i,
            "venue": "Mock Journal",
            "type": "journal-article",
        }
        for i in (1, 2)
    ]


def _parse_sources(raw: str | None) -> tuple[list[str], list[str]]:
    names = [normalize_source(s) for s in (raw or "").replace(";", ",").split(",") if s.strip()]
    names = names or list(DEFAULT_SOURCES)
    seen: list[str] = []
    for n in names:
        if n not in seen:
            seen.append(n)
    bad = [n for n in seen if n not in SUPPORTED_SOURCES]
    return seen, bad


def main() -> int:
    p = argparse.ArgumentParser(description="Search external sources and save hits to the library")
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    p.add_argument("--query", required=True)
    p.add_argument("--sources", default="", help="Comma-separated; default pubmed,europepmc")
    p.add_argument("--date_from", default="")
    p.add_argument("--date_to", default="")
    p.add_argument("--limit", type=int, default=20, help="Hits per source (1-50)")
    p.add_argument("--run_id", default="", help="Optional id; a new one is made when empty")
    p.add_argument("--mock", action="store_true", help="Offline fake hits for tests")
    p.add_argument("--skip_ingest", action="store_true", help="Write hits.jsonl only")
    args = p.parse_args()

    sources, bad = _parse_sources(args.sources)
    if bad:
        return fail(f"unsupported sources: {bad}", supported=list(SUPPORTED_SOURCES))
    if len(sources) > MAX_SOURCES:
        return fail(f"at most {MAX_SOURCES} sources per search")
    limit = max(1, min(int(args.limit), MAX_LIMIT))

    conn, err = open_project_db(args.project_root)
    if err:
        return fail(err)
    try:
        col, err = resolve_collection(conn, args.collection_id or None)
        if err:
            return fail(err)
        collection_id, col_path = col["id"], Path(col["path"])
    finally:
        conn.close()

    run_id = args.run_id.strip() or f"run-{uuid.uuid4().hex[:10]}"
    run_dir = col_path / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    protocol = build_frozen_protocol(
        run_id=run_id,
        collection_id=collection_id,
        query=args.query,
        sources=sources,
        date_from=args.date_from or None,
        date_to=args.date_to or None,
        limit=limit,
    )

    all_hits: list[dict] = []
    per_source: list[dict] = []
    executed_at = utc_now()
    for src in sources:
        if args.mock:
            hits, meta = _mock_hits(args.query, src), {"provider": "mock", "source": src}
        else:
            hits, meta = search_source(
                src, args.query, limit=limit,
                date_from=args.date_from or None, date_to=args.date_to or None,
            )
        failed = meta.get("provider") == "error"
        entry = {"source": src, "n": len(hits), "ok": not failed}
        if failed:
            entry["error"] = meta.get("message")
        per_source.append(entry)
        all_hits.extend(hits)

    hits_path = run_dir / "hits.jsonl"
    with hits_path.open("w", encoding="utf-8") as f:
        for h in all_hits:
            f.write(json.dumps(h, ensure_ascii=False) + "\n")
    protocol.update(
        executed_at=executed_at,
        per_source=per_source,
        n_hits=len(all_hits),
        status=RUN_STATUS_SEARCHED,
        hits_path=str(hits_path),
    )
    (run_dir / "run.json").write_text(
        json.dumps(protocol, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    db_path = Path(args.project_root).resolve() / "library.sqlite"
    conn = open_db(db_path)
    try:
        conn.execute(
            """
            INSERT OR REPLACE INTO runs(
              id, collection_id, query, sources_json, config_json,
              started_at, finished_at, hits_path, status, coverage_path, gate_json
            ) VALUES (?, ?, ?, ?, ?, ?, NULL, NULL, ?, NULL, NULL)
            """,
            (
                run_id,
                collection_id,
                args.query,
                json.dumps(sources),
                json.dumps(protocol, ensure_ascii=False),
                executed_at,
                RUN_STATUS_SEARCHED,
            ),
        )
        set_run_status(conn, run_id, RUN_STATUS_SEARCHED, hits_path=str(hits_path))
        conn.commit()
    finally:
        conn.close()

    ingest = None
    if not args.skip_ingest and all_hits:
        ingest = ingest_hits_file(
            db_path=db_path, collection_id=collection_id, run_id=run_id, hits_path=hits_path
        )

    failed_sources = [s["source"] for s in per_source if not s["ok"]]
    emit(
        {
            "ok": True,
            "run_id": run_id,
            "collection_id": collection_id,
            "query": args.query,
            "mode": "mock" if args.mock else ("partial" if failed_sources else "live"),
            "n_hits": len(all_hits),
            "per_source": per_source,
            "hits_path": str(hits_path),
            "ingest": ingest,
            "warnings": [f"source failed: {s}" for s in failed_sources],
        }
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
