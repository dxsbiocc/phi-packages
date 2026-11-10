#!/usr/bin/env python3
"""Offline enrich smoke (--mock Zheng fixture)."""
from __future__ import annotations

import json
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

from db import init_db, open_db  # noqa: E402
from enrich import enrich_publications, needs_enrich  # noqa: E402
from ingest import ingest_records, load_records  # noqa: E402

ENRICH = ROOT / "skills" / "manage-library" / "scripts" / "lib" / "enrich.py"
INGEST = ROOT / "skills" / "manage-library" / "scripts" / "ingest.py"
COL = ROOT / "skills" / "manage-library" / "scripts" / "collection.py"
ZHENG_TITLE = "Landscape of Infiltrating T Cells in Liver Cancer Revealed by Single-Cell Sequencing"


def _make_project(tmp: Path) -> tuple[Path, str]:
    db = tmp / "library.sqlite"
    init_db(db)
    r = subprocess.run(
        [
            sys.executable,
            str(COL),
            "--project_root",
            str(tmp),
            "--action",
            "create",
            "--name",
            "enrich-test",
            "--no_wiki",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    out = json.loads(r.stdout.strip().splitlines()[-1])
    return db, out["collection_id"]


class NeedsEnrichTests(unittest.TestCase):
    def test_shell_needs(self):
        self.assertTrue(
            needs_enrich({"title": "10.1016/j.cell.2017.05.035", "doi": "10.1016/j.cell.2017.05.035", "abstract": ""})
        )
        self.assertFalse(
            needs_enrich({"title": "Real Title", "abstract": "x" * 50, "doi": "10.1/x"})
        )


class EnrichMockTests(unittest.TestCase):
    def test_enrich_mock_standalone(self):
        with tempfile.TemporaryDirectory() as td:
            proj = Path(td)
            db, cid = _make_project(proj)
            ids = proj / "zheng.ids"
            ids.write_text("doi:10.1016/j.cell.2017.05.035\npmid:28622514\n", encoding="utf-8")
            records = load_records(ids, fmt="ids")
            ing = ingest_records(db_path=db, collection_id=cid, records=records)
            self.assertTrue(ing["ok"])
            # still a shell
            conn = open_db(db)
            row = dict(conn.execute("SELECT * FROM publications").fetchone())
            conn.close()
            self.assertTrue(needs_enrich(row))
            en = enrich_publications(db_path=db, collection_id=cid, mock=True)
            self.assertTrue(en["ok"])
            self.assertGreaterEqual(en["n_enriched"], 1)
            conn = open_db(db)
            row = dict(conn.execute("SELECT * FROM publications").fetchone())
            conn.close()
            self.assertIn("Landscape of Infiltrating T Cells", row["title"] or "")
            self.assertGreater(len(row["abstract"] or ""), 200)
            self.assertEqual(row["pmid"], "28622514")
            self.assertEqual((row["doi"] or "").lower(), "10.1016/j.cell.2017.05.035")

    def test_cli_mock(self):
        with tempfile.TemporaryDirectory() as td:
            proj = Path(td)
            db, cid = _make_project(proj)
            ids = proj / "zheng.ids"
            ids.write_text("doi:10.1016/j.cell.2017.05.035\npmid:28622514\n", encoding="utf-8")
            r = subprocess.run(
                [
                    sys.executable,
                    str(INGEST),
                    "--project_root",
                    str(proj),
                    "--collection_id",
                    cid,
                    "--input_path",
                    str(ids),
                    "--format",
                    "ids",
                    "--enrich_mock",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            out = json.loads(r.stdout.strip().splitlines()[-1])
            self.assertTrue(out["ok"])
            self.assertGreaterEqual(out.get("enrich", {}).get("n_enriched", 0), 1)
            conn = sqlite3.connect(db)
            title, abstract = conn.execute("SELECT title, abstract FROM publications").fetchone()
            conn.close()
            self.assertIn("Landscape", title)
            self.assertGreater(len(abstract or ""), 200)


if __name__ == "__main__":
    unittest.main()
