#!/usr/bin/env python3
"""Unit tests for normalize_* and resolve_or_create (no network)."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

from db import (  # noqa: E402
    SCHEMA_VERSION,
    apply_id_merges,
    init_db,
    normalize_arxiv,
    normalize_doi,
    normalize_openalex,
    normalize_pmid,
    normalize_pmcid,
    open_db,
    plan_id_merges,
    resolve_or_create,
    soft_fingerprint,
)


class NormalizeTests(unittest.TestCase):
    def test_doi(self):
        self.assertEqual(normalize_doi("https://doi.org/10.1016/j.cell.2017.05.035"), "10.1016/j.cell.2017.05.035")
        self.assertEqual(normalize_doi("DOI:10.1016/J.CELL.2017.05.035"), "10.1016/j.cell.2017.05.035")
        self.assertIsNone(normalize_doi(""))
        self.assertIsNone(normalize_doi("not-a-doi"))

    def test_pmid(self):
        self.assertEqual(normalize_pmid("PMID:28622514"), "28622514")
        self.assertEqual(normalize_pmid("28622514"), "28622514")

    def test_pmcid(self):
        self.assertEqual(normalize_pmcid("pmc12345"), "PMC12345")
        self.assertEqual(normalize_pmcid("PMCID:PMC12345"), "PMC12345")

    def test_arxiv(self):
        self.assertEqual(normalize_arxiv("https://arxiv.org/abs/1234.5678"), "1234.5678")
        self.assertEqual(normalize_arxiv("arXiv:1234.5678"), "1234.5678")

    def test_openalex(self):
        self.assertEqual(normalize_openalex("https://openalex.org/W123"), "W123")
        self.assertEqual(normalize_openalex("openalex:W123"), "W123")


class ResolveTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db = Path(self.tmp.name) / "library.sqlite"
        init_db(self.db)
        self.conn = open_db(self.db)

    def tearDown(self):
        self.conn.close()
        self.tmp.cleanup()

    def test_insert_then_merge_by_doi(self):
        r1 = resolve_or_create(
            self.conn,
            {
                "doi": "10.1000/test.1",
                "title": "First Title",
                "authors": ["A, B"],
                "year": 2020,
            },
        )
        self.conn.commit()
        self.assertEqual(r1["action"], "inserted")
        r2 = resolve_or_create(
            self.conn,
            {
                "doi": "https://doi.org/10.1000/TEST.1",
                "pmid": "111",
                "title": "First Title Updated Abstract Carrier",
                "abstract": "Hello",
                "authors": ["A, B"],
                "year": 2020,
            },
        )
        self.conn.commit()
        self.assertEqual(r2["action"], "merged")
        self.assertEqual(r2["publication_id"], r1["publication_id"])
        row = self.conn.execute(
            "SELECT pmid, abstract FROM publications WHERE id = ?", (r1["publication_id"],)
        ).fetchone()
        self.assertEqual(row["pmid"], "111")
        self.assertEqual(row["abstract"], "Hello")

    def test_soft_conflict_needs_review(self):
        soft = soft_fingerprint("Soft Fingerprint Collision Example About Exhaustion", ["Doe, Jane"], 2019)
        self.assertIsNotNone(soft)
        r1 = resolve_or_create(
            self.conn,
            {
                "doi": "10.1000/soft.a",
                "title": "Soft Fingerprint Collision Example About Exhaustion",
                "authors": ["Doe, Jane"],
                "year": 2019,
            },
        )
        self.conn.commit()
        self.assertEqual(r1["action"], "inserted")
        r2 = resolve_or_create(
            self.conn,
            {
                "doi": "10.1000/soft.b",
                "title": "Soft Fingerprint Collision Example About Exhaustion",
                "authors": ["Doe, Jane"],
                "year": 2019,
            },
        )
        self.conn.commit()
        # Different strong DOI but same soft fingerprint → needs_review, no merge
        self.assertEqual(r2["action"], "needs_review")
        self.assertNotEqual(r2["publication_id"], r1["publication_id"])
        row = self.conn.execute(
            "SELECT needs_review FROM publications WHERE id = ?", (r2["publication_id"],)
        ).fetchone()
        self.assertEqual(row["needs_review"], 1)

    def test_conflicting_strong_ids(self):
        resolve_or_create(self.conn, {"doi": "10.1000/a", "title": "A", "year": 2021, "authors": ["X"]})
        resolve_or_create(self.conn, {"pmid": "222", "title": "B", "year": 2022, "authors": ["Y"]})
        self.conn.commit()
        # Draft claims both ids that already map to different pubs
        r = resolve_or_create(
            self.conn,
            {"doi": "10.1000/a", "pmid": "222", "title": "Conflict", "year": 2023, "authors": ["Z"]},
        )
        self.conn.commit()
        self.assertEqual(r["action"], "needs_review")
        self.assertEqual(r.get("reason"), "conflicting_strong_ids")

    def test_complementary_doi_pmid_shells_auto_merge(self):
        """doi-only shell + pmid-only shell + joint hit → one doi: canonical."""
        r1 = resolve_or_create(
            self.conn,
            {"doi": "10.1016/j.cell.2017.05.035", "title": "10.1016/j.cell.2017.05.035"},
        )
        self.conn.commit()
        self.assertEqual(r1["action"], "inserted")
        r2 = resolve_or_create(
            self.conn,
            {"pmid": "28622514", "title": "28622514"},
        )
        self.conn.commit()
        self.assertEqual(r2["action"], "inserted")
        self.assertNotEqual(r1["publication_id"], r2["publication_id"])
        r3 = resolve_or_create(
            self.conn,
            {
                "doi": "10.1016/j.cell.2017.05.035",
                "pmid": "28622514",
                "title": "Landscape of Infiltrating T Cells in Liver Cancer",
                "abstract": "LAYN and exhausted CD8 T cells in HCC.",
                "year": 2017,
                "authors": ["Zheng, C"],
            },
        )
        self.conn.commit()
        self.assertEqual(r3["action"], "merged")
        self.assertEqual(r3["publication_id"], "doi:10.1016/j.cell.2017.05.035")
        self.assertEqual(r3.get("reason"), "complementary_strong_ids")
        n = self.conn.execute("SELECT COUNT(*) c FROM publications").fetchone()["c"]
        self.assertEqual(n, 1)
        row = self.conn.execute(
            "SELECT doi, pmid, needs_review, abstract, title FROM publications WHERE id=?",
            ("doi:10.1016/j.cell.2017.05.035",),
        ).fetchone()
        self.assertEqual(row["pmid"], "28622514")
        self.assertEqual(row["needs_review"], 0)
        self.assertIn("LAYN", row["abstract"])
        self.assertIn("Landscape", row["title"])

    def test_repair_plan_and_apply_on_conflict_cluster(self):
        """Existing shell+conflict rows are repaired by apply_id_merges."""
        # Simulate pre-0.2.2 state: two shells + conflict row (manual insert path)
        resolve_or_create(self.conn, {"doi": "10.1000/repair.1", "title": "10.1000/repair.1"})
        resolve_or_create(self.conn, {"pmid": "999001", "title": "999001"})
        self.conn.commit()
        # Force old conflict behavior by inserting conflict-like third row sharing ids
        # Use resolve with both — should auto-merge now; so build cluster manually for repair test
        # Clear and build legacy shape:
        for table in ("aliases", "pubs_fts", "publications"):
            try:
                self.conn.execute(f"DELETE FROM {table}")
            except Exception:
                pass
        self.conn.commit()
        from db import _insert_pub, utc_now
        now = utc_now()
        _insert_pub(
            self.conn,
            "doi:10.1000/repair.1",
            {"doi": "10.1000/repair.1", "pmid": None, "pmcid": None, "arxiv": None, "openalex": None},
            {"title": "10.1000/repair.1"},
            None,
            now,
            needs_review=0,
        )
        self.conn.execute(
            "INSERT INTO aliases(alias, publication_id) VALUES (?, ?)",
            ("doi:10.1000/repair.1", "doi:10.1000/repair.1"),
        )
        _insert_pub(
            self.conn,
            "pmid:999001",
            {"doi": None, "pmid": "999001", "pmcid": None, "arxiv": None, "openalex": None},
            {"title": "999001"},
            None,
            now,
            needs_review=0,
        )
        self.conn.execute(
            "INSERT INTO aliases(alias, publication_id) VALUES (?, ?)",
            ("pmid:999001", "pmid:999001"),
        )
        conflict = "doi:10.1000/repair.1#conflict-deadbeef"
        _insert_pub(
            self.conn,
            conflict,
            {"doi": "10.1000/repair.1", "pmid": "999001", "pmcid": None, "arxiv": None, "openalex": None},
            {
                "title": "Repair Target Paper About Exhaustion",
                "abstract": "Full abstract here.",
                "year": 2020,
                "authors": ["A"],
            },
            None,
            now,
            needs_review=1,
        )
        self.conn.commit()
        plan = plan_id_merges(self.conn, only_high_confidence=True)
        self.assertGreaterEqual(plan["n_actionable"], 1)
        out = apply_id_merges(self.conn, only_high_confidence=True)
        self.conn.commit()
        self.assertGreaterEqual(out["n_applied"], 1)
        n = self.conn.execute("SELECT COUNT(*) c FROM publications").fetchone()["c"]
        self.assertEqual(n, 1)
        row = self.conn.execute("SELECT * FROM publications").fetchone()
        self.assertEqual(row["id"], "doi:10.1000/repair.1")
        self.assertEqual(row["pmid"], "999001")
        self.assertEqual(row["needs_review"], 0)
        self.assertIn("Repair Target", row["title"])


    def test_pmid_then_doi_upgrades_pk(self):
        """pmid-only insert then doi+pmid → single row with doi:… PK."""
        r1 = resolve_or_create(
            self.conn,
            {
                "pmid": "28622514",
                "title": "Landscape of Infiltrating T Cells in Liver Cancer",
                "authors": ["Zheng, C"],
                "year": 2017,
            },
        )
        self.conn.commit()
        self.assertEqual(r1["action"], "inserted")
        self.assertEqual(r1["publication_id"], "pmid:28622514")
        r2 = resolve_or_create(
            self.conn,
            {
                "doi": "10.1016/j.cell.2017.05.035",
                "pmid": "28622514",
                "title": "Landscape of Infiltrating T Cells in Liver Cancer",
                "abstract": "LAYN and exhausted CD8 T cells.",
                "authors": ["Zheng, C"],
                "year": 2017,
            },
        )
        self.conn.commit()
        self.assertEqual(r2["action"], "merged")
        self.assertEqual(r2["publication_id"], "doi:10.1016/j.cell.2017.05.035")
        self.assertEqual(r2.get("previous_id"), "pmid:28622514")
        self.assertEqual(r2.get("reason"), "canonical_id_upgrade")
        n = self.conn.execute("SELECT COUNT(*) c FROM publications").fetchone()["c"]
        self.assertEqual(n, 1)
        row = self.conn.execute("SELECT * FROM publications").fetchone()
        self.assertEqual(row["id"], "doi:10.1016/j.cell.2017.05.035")
        self.assertEqual(row["pmid"], "28622514")
        self.assertEqual(row["doi"], "10.1016/j.cell.2017.05.035")
        aliases = {
            r["alias"]: r["publication_id"]
            for r in self.conn.execute("SELECT alias, publication_id FROM aliases")
        }
        self.assertEqual(aliases.get("doi:10.1016/j.cell.2017.05.035"), row["id"])
        self.assertEqual(aliases.get("pmid:28622514"), row["id"])

    def test_pmcid_then_pmid_upgrades_pk(self):
        """pmcid-only then pmid+pmcid → PK upgrades to pmid:…"""
        r1 = resolve_or_create(
            self.conn,
            {"pmcid": "PMC5506419", "title": "PMC shell only", "year": 2017, "authors": ["Z"]},
        )
        self.conn.commit()
        self.assertEqual(r1["publication_id"], "pmcid:PMC5506419")
        r2 = resolve_or_create(
            self.conn,
            {
                "pmid": "28622514",
                "pmcid": "PMC5506419",
                "title": "Real title after PMID arrives",
                "year": 2017,
                "authors": ["Zheng, C"],
            },
        )
        self.conn.commit()
        self.assertEqual(r2["action"], "merged")
        self.assertEqual(r2["publication_id"], "pmid:28622514")
        self.assertEqual(r2.get("previous_id"), "pmcid:PMC5506419")
        n = self.conn.execute("SELECT COUNT(*) c FROM publications").fetchone()["c"]
        self.assertEqual(n, 1)
        row = self.conn.execute("SELECT pmid, pmcid FROM publications WHERE id=?", ("pmid:28622514",)).fetchone()
        self.assertEqual(row["pmid"], "28622514")
        self.assertEqual(row["pmcid"], "PMC5506419")
        aliases = {
            r["alias"]: r["publication_id"]
            for r in self.conn.execute("SELECT alias, publication_id FROM aliases")
        }
        self.assertEqual(aliases.get("pmid:28622514"), "pmid:28622514")
        self.assertEqual(aliases.get("pmcid:PMC5506419"), "pmid:28622514")

    def test_arxiv_then_doi_upgrades_pk(self):
        r1 = resolve_or_create(
            self.conn,
            {"arxiv": "1706.01814", "title": "Arxiv preprint", "year": 2017, "authors": ["A"]},
        )
        self.conn.commit()
        self.assertEqual(r1["publication_id"], "arxiv:1706.01814")
        r2 = resolve_or_create(
            self.conn,
            {
                "doi": "10.1000/arxiv.upgrade.1",
                "arxiv": "1706.01814",
                "title": "Published version",
                "year": 2018,
                "authors": ["A"],
            },
        )
        self.conn.commit()
        self.assertEqual(r2["action"], "merged")
        self.assertEqual(r2["publication_id"], "doi:10.1000/arxiv.upgrade.1")
        self.assertEqual(r2.get("previous_id"), "arxiv:1706.01814")
        n = self.conn.execute("SELECT COUNT(*) c FROM publications").fetchone()["c"]
        self.assertEqual(n, 1)

    def test_schema_unique_doi_index(self):
        ver = self.conn.execute(
            "SELECT value FROM meta WHERE key='schema_version'"
        ).fetchone()["value"]
        self.assertEqual(ver, str(SCHEMA_VERSION))
        self.assertGreaterEqual(SCHEMA_VERSION, 2)
        resolve_or_create(
            self.conn,
            {"doi": "10.1000/unique.gate", "title": "One", "year": 2020, "authors": ["A"]},
        )
        self.conn.commit()
        # resolve path must not create a second row for the same DOI
        r2 = resolve_or_create(
            self.conn,
            {"doi": "10.1000/UNIQUE.GATE", "title": "One again", "year": 2020, "authors": ["A"]},
        )
        self.conn.commit()
        self.assertEqual(r2["action"], "merged")
        n = self.conn.execute("SELECT COUNT(*) c FROM publications").fetchone()["c"]
        self.assertEqual(n, 1)
        # Raw SQL duplicate DOI should fail under UNIQUE index
        import sqlite3
        with self.assertRaises(sqlite3.IntegrityError):
            self.conn.execute(
                "INSERT INTO publications(id, doi, title, needs_review, created_at, updated_at) "
                "VALUES ('doi:dup-raw', '10.1000/unique.gate', 'dup', 0, 't', 't')"
            )
            self.conn.commit()
        self.conn.rollback()




if __name__ == "__main__":
    unittest.main()
