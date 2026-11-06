#!/usr/bin/env python3
"""Write runs/<run_id>/run.yaml + run.json with frozen search protocol (0.5.2)."""
from __future__ import annotations

import argparse
import json
import sys
import uuid
from pathlib import Path

_LIB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
sys.path.insert(0, str(_LIB))
from db import (  # noqa: E402
    RUN_STATUS_CONFIGURED,
    build_frozen_protocol,
    open_db,
    utc_now,
)

ALLOWED_SOURCES = [
    "pubmed",
    "europepmc",
    "openalex",
    "crossref",
    "semantic_scholar",
    "arxiv",
    "biorxiv",
    "medrxiv",
    "pmc",
    "doaj",
    "datacite",
    "zenodo",
]


def _load_id_list(raw: str) -> list[str]:
    """Parse comma list or @file (one ID per line / comma-separated)."""
    s = (raw or "").strip()
    if not s:
        return []
    if s.startswith("@"):
        path = Path(s[1:]).expanduser()
        body = path.read_text(encoding="utf-8")
        items: list[str] = []
        for line in body.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            items.extend(x.strip() for x in line.split(",") if x.strip())
        return items
    return [x.strip() for x in s.split(",") if x.strip()]



def main() -> int:
    p = argparse.ArgumentParser(description="Freeze a literature search run protocol")
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    p.add_argument("--query", required=True)
    p.add_argument(
        "--sources",
        default="pubmed,europepmc",
        help="Comma-separated open sources (draft allows up to 8)",
    )
    p.add_argument("--date_from", default="")
    p.add_argument("--date_to", default="")
    p.add_argument("--type", default="any", dest="pub_type")
    p.add_argument(
        "--limit",
        type=int,
        default=20,
        help="Per-source hit limit (frozen as limit + limit_per_source)",
    )
    p.add_argument("--run_id", default="")
    p.add_argument(
        "--gold_set",
        default="",
        help="Comma-separated verified strong IDs (DOI preferred; pmid/arxiv ok). "
        "Guidance 5–10; <5 marks evidence_weak (not hard-fail). Or @path to a text file.",
    )
    p.add_argument(
        "--out_of_window_controls",
        default="",
        help="Comma-separated 1–2 real papers outside the date window (prove date filter). Or @path.",
    )
    args = p.parse_args()

    alias = {
        "s2": "semantic_scholar",
        "semantic-scholar": "semantic_scholar",
        "epmc": "europepmc",
        "europe_pmc": "europepmc",
        "pubmed_central": "pmc",
        "pubmed-central": "pmc",
        "data_cite": "datacite",
        "data-cite": "datacite",
    }
    sources = [
        alias.get(s.strip().lower(), s.strip().lower())
        for s in args.sources.split(",")
        if s.strip()
    ]
    bad = [s for s in sources if s not in ALLOWED_SOURCES]
    if bad:
        print(
            json.dumps(
                {"ok": False, "error": f"unknown sources: {bad}", "allowed": ALLOWED_SOURCES}
            )
        )
        return 1
    if len(sources) > 8:
        print(json.dumps({"ok": False, "error": "at most 8 sources per run"}))
        return 1
    if not sources:
        sources = ["pubmed"]

    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    if not db_path.exists():
        print(json.dumps({"ok": False, "error": f"missing {db_path}; create collection first"}))
        return 1

    conn = open_db(db_path)
    try:
        col = conn.execute(
            "SELECT id, path FROM collections WHERE id = ?", (args.collection_id,)
        ).fetchone()
        if not col:
            print(
                json.dumps(
                    {"ok": False, "error": f"unknown collection_id {args.collection_id}"}
                )
            )
            return 1
        col_path = Path(col["path"])
    finally:
        conn.close()

    run_id = args.run_id.strip() or f"run-{uuid.uuid4().hex[:10]}"
    run_dir = col_path / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    limit = max(1, min(args.limit, 50))
    gold_set = _load_id_list(args.gold_set)
    oow = _load_id_list(args.out_of_window_controls)
    config = build_frozen_protocol(
        run_id=run_id,
        collection_id=args.collection_id,
        query=args.query,
        sources=sources,
        date_from=args.date_from or None,
        date_to=args.date_to or None,
        pub_type=args.pub_type,
        limit=limit,
        gold_set=gold_set,
        out_of_window_controls=oow,
    )
    config["status"] = RUN_STATUS_CONFIGURED

    lines = [
        f"protocol_version: {config['protocol_version']}",
        f"run_id: {config['run_id']}",
        f"collection_id: {config['collection_id']}",
        f"query: {json.dumps(config['query'], ensure_ascii=False)}",
        "sources:",
    ]
    for s in sources:
        lines.append(f"  - {s}")
    lines.append(f"date_from: {config['date_from'] if config['date_from'] else 'null'}")
    lines.append(f"date_to: {config['date_to'] if config['date_to'] else 'null'}")
    lines.append(f"type: {config['type']}")
    lines.append(f"limit: {config['limit']}")
    lines.append(f"limit_per_source: {config['limit_per_source']}")
    lines.append("provider_notes:")
    for s, note in (config.get("provider_notes") or {}).items():
        lines.append(f"  {s}: {json.dumps(note, ensure_ascii=False)}")
    lines.append("gold_set:")
    for g in config.get("gold_set") or []:
        lines.append(f"  - {g}")
    if not (config.get("gold_set") or []):
        lines.append("  []")
    lines.append(f"evidence_weak_if_gold_lt: {config.get('gold_set_min_guidance', 5)}")
    lines.append("out_of_window_controls:")
    for g in config.get("out_of_window_controls") or []:
        lines.append(f"  - {g}")
    if not (config.get("out_of_window_controls") or []):
        lines.append("  []")
    lines.append(f"frozen_at: {config['frozen_at']}")
    lines.append(f"created_at: {config['created_at']}")
    lines.append(f"status: {config['status']}")
    lines.append("no_commercial_sources: true")
    yaml_path = run_dir / "run.yaml"
    yaml_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (run_dir / "run.json").write_text(
        json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    conn = open_db(db_path)
    try:
        conn.execute(
            """
            INSERT OR REPLACE INTO runs(
              id, collection_id, query, sources_json, config_json,
              started_at, finished_at, hits_path, status, coverage_path, gate_json
            )
            VALUES (?, ?, ?, ?, ?, ?, NULL, NULL, ?, NULL, NULL)
            """,
            (
                run_id,
                args.collection_id,
                args.query,
                json.dumps(sources),
                json.dumps(config, ensure_ascii=False),
                utc_now(),
                RUN_STATUS_CONFIGURED,
            ),
        )
        conn.commit()
    finally:
        conn.close()

    print(
        json.dumps(
            {
                "ok": True,
                "mode": "live",
                "run_id": run_id,
                "run_yaml": str(yaml_path),
                "run_json": str(run_dir / "run.json"),
                "config": config,
                "protocol_frozen": True,
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
