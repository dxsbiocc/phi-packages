#!/usr/bin/env python3
"""Smoke helper: DOI shell + PMID shell + joint hit → one canonical pub."""
from __future__ import annotations

import argparse
import json
import sqlite3
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    args = p.parse_args()
    demo = Path(args.project_root).resolve()
    cid = args.collection_id
    py = sys.executable

    ids = demo / "zheng_ids.txt"
    ids.write_text("doi:10.1016/j.cell.2017.05.035\nPMID:28622514\n", encoding="utf-8")
    subprocess.check_call(
        [
            py,
            str(ROOT / "skills/literature-library/scripts/ingest.py"),
            "--project_root",
            str(demo),
            "--collection_id",
            cid,
            "--input_path",
            str(ids),
            "--format",
            "ids",
        ]
    )
    hits = demo / "zheng_joint.jsonl"
    hits.write_text(
        json.dumps(
            {
                "doi": "10.1016/j.cell.2017.05.035",
                "pmid": "28622514",
                "title": "Landscape of Infiltrating T Cells in Liver Cancer Revealed by Single-Cell Sequencing",
                "abstract": (
                    "Here we report single-cell RNA sequencing of T cells from HCC. "
                    "Exhausted CD8+ T cells and LAYN."
                ),
                "year": 2017,
                "authors": ["Zheng, C"],
                "source": "test",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    out = subprocess.check_output(
        [
            py,
            str(ROOT / "skills/literature-library/scripts/ingest.py"),
            "--project_root",
            str(demo),
            "--collection_id",
            cid,
            "--input_path",
            str(hits),
            "--format",
            "hits",
        ],
        text=True,
    )
    o = json.loads(out)
    assert o["ok"], o

    db = sqlite3.connect(str(demo / "library.sqlite"))
    rows = db.execute(
        """
        SELECT id, doi, pmid, needs_review, length(ifnull(abstract,''))
        FROM publications
        WHERE doi='10.1016/j.cell.2017.05.035'
           OR pmid='28622514'
           OR id LIKE '%10.1016/j.cell.2017.05.035%'
        """
    ).fetchall()
    print("zheng rows", rows)
    assert len(rows) == 1, rows
    assert rows[0][0] == "doi:10.1016/j.cell.2017.05.035", rows
    assert rows[0][2] == "28622514" and rows[0][3] == 0 and rows[0][4] > 20, rows

    plan = subprocess.check_output(
        [py, str(ROOT / "skills/literature-library/scripts/merge_ids_cli.py"), "--project_root", str(demo)],
        text=True,
    )
    pobj = json.loads(plan)
    assert pobj["ok"] and pobj["mode"] == "dry-run", pobj
    print("merge dry-run actionable", pobj.get("n_actionable"))
    db.close()
    print("MERGE SMOKE OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
