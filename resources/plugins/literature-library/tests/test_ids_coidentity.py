#!/usr/bin/env python3
"""Ids-file co-identity: dual DOI+PMID → one publication; unrelated stay separate."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "lit-ingest" / "scripts"))
sys.path.insert(0, str(ROOT / "skills" / "lit-db" / "scripts"))
sys.path.insert(0, str(ROOT / "skills" / "lit-collection" / "scripts"))

from ingest import _parse_id_list, ingest_records  # noqa: E402
from db import init_db, open_db  # noqa: E402


class ParseCoIdentityTests(unittest.TestCase):
    def test_zheng_two_lines(self):
        rows = _parse_id_list("doi:10.1016/j.cell.2017.05.035\npmid:28622514\n")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["doi"], "10.1016/j.cell.2017.05.035")
        self.assertEqual(rows[0]["pmid"], "28622514")
        self.assertEqual(rows[0].get("co_identity"), "ids_file")

    def test_same_line_and_tab(self):
        a = _parse_id_list("doi:10.1016/j.cell.2017.05.035 pmid:28622514\n")
        b = _parse_id_list("10.1016/j.cell.2017.05.035\t28622514\n")
        self.assertEqual(len(a), 1)
        self.assertEqual(len(b), 1)
        self.assertEqual(a[0]["pmid"], "28622514")
        self.assertEqual(b[0]["pmid"], "28622514")

    def test_alternating_two_papers(self):
        rows = _parse_id_list(
            "doi:10.1000/a\npmid:111\ndoi:10.1000/b\npmid:222\n"
        )
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["doi"], "10.1000/a")
        self.assertEqual(rows[0]["pmid"], "111")
        self.assertEqual(rows[1]["doi"], "10.1000/b")
        self.assertEqual(rows[1]["pmid"], "222")

    def test_blank_groups(self):
        rows = _parse_id_list("doi:10.1000/a\npmid:111\n\ndoi:10.1000/b\npmid:222\n")
        self.assertEqual(len(rows), 2)

    def test_block_layout_not_cross_merged(self):
        rows = _parse_id_list(
            "doi:10.1000/a\ndoi:10.1000/b\npmid:111\npmid:222\n"
        )
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(r.get("pmid") is None or r.get("doi") is None for r in rows))


class IngestCoIdentityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.db = self.root / "library.sqlite"
        init_db(self.db)
        # Minimal collection row matching schema
        conn = open_db(self.db)
        cols = [r[1] for r in conn.execute("PRAGMA table_info(collections)").fetchall()]
        path = str(self.root / "collections" / "col-t")
        Path(path).mkdir(parents=True, exist_ok=True)
        fields = {"id": "col-t", "name": "t", "path": path}
        # fill required columns with sensible defaults
        values = []
        names = []
        for c in cols:
            names.append(c)
            if c in fields:
                values.append(fields[c])
            elif c.endswith("_at"):
                values.append("2026-01-01T00:00:00Z")
            else:
                values.append(None)
        placeholders = ",".join("?" for _ in names)
        conn.execute(
            f"INSERT INTO collections ({','.join(names)}) VALUES ({placeholders})",
            values,
        )
        conn.commit()
        conn.close()

    def tearDown(self):
        self.tmp.cleanup()

    def test_zheng_ids_only_one_pub(self):
        records = _parse_id_list(
            "doi:10.1016/j.cell.2017.05.035\npmid:28622514\n"
        )
        self.assertEqual(len(records), 1)
        out = ingest_records(
            db_path=self.db, collection_id="col-t", records=records
        )
        self.assertTrue(out["ok"], out)
        conn = open_db(self.db)
        rows = conn.execute(
            "SELECT id, doi, pmid FROM publications ORDER BY id"
        ).fetchall()
        conn.close()
        self.assertEqual(len(rows), 1, rows)
        self.assertEqual(rows[0][0], "doi:10.1016/j.cell.2017.05.035")
        self.assertEqual(rows[0][1], "10.1016/j.cell.2017.05.035")
        self.assertEqual(rows[0][2], "28622514")


    def test_ids_then_bib_enriches_title(self):
        records = _parse_id_list(
            "doi:10.1016/j.cell.2017.05.035\npmid:28622514\n"
        )
        ingest_records(db_path=self.db, collection_id="col-t", records=records)
        rich = [{
            "doi": "10.1016/j.cell.2017.05.035",
            "pmid": "28622514",
            "title": "Landscape of Infiltrating T Cells in Liver Cancer Revealed by Single-Cell Sequencing",
            "abstract": "Deep single-cell RNA sequencing on T cells from HCC patients.",
            "year": 2017,
            "source": "bib",
        }]
        out = ingest_records(db_path=self.db, collection_id="col-t", records=rich)
        self.assertTrue(out["ok"], out)
        conn = open_db(self.db)
        row = conn.execute(
            "SELECT id, title, pmid, length(abstract) FROM publications"
        ).fetchone()
        conn.close()
        self.assertEqual(row[0], "doi:10.1016/j.cell.2017.05.035")
        self.assertTrue(str(row[1]).startswith("Landscape"), row)
        self.assertEqual(row[2], "28622514")
        self.assertGreater(row[3], 20)

    def test_two_papers_stay_separate(self):
        records = _parse_id_list(
            "doi:10.1000/paper.a\npmid:1001\n\ndoi:10.1000/paper.b\npmid:1002\n"
        )
        out = ingest_records(
            db_path=self.db, collection_id="col-t", records=records
        )
        self.assertTrue(out["ok"], out)
        conn = open_db(self.db)
        n = conn.execute("SELECT COUNT(*) FROM publications").fetchone()[0]
        conn.close()
        self.assertEqual(n, 2)


if __name__ == "__main__":
    unittest.main()
