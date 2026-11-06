#!/usr/bin/env python3
"""Offline tests for human_queue, coverage_report, publication_versions + schema v4 gate cols."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "literature-library" / "scripts"))

from db import (  # noqa: E402
    SCHEMA_VERSION,
    apply_version_detection,
    confirm_version_link,
    coverage_report,
    detect_version_candidates,
    init_db,
    link_to_collection,
    list_human_queue,
    merge_preprint_journal,
    open_db,
    propose_version_link,
    record_run_hit,
    resolve_or_create,
    utc_now,
)


class SchemaV3Tests(unittest.TestCase):
    def test_schema_version_4_and_versions_table(self):
        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "library.sqlite"
            info = init_db(db)
            self.assertGreaterEqual(info["schema_version"], 4)
            self.assertGreaterEqual(SCHEMA_VERSION, 4)
            conn = open_db(db)
            try:
                row = conn.execute(
                    "SELECT value FROM meta WHERE key='schema_version'"
                ).fetchone()
                self.assertEqual(row["value"], str(SCHEMA_VERSION))
                cols = {r[1] for r in conn.execute("PRAGMA table_info(publications)")}
                self.assertIn("review_reasons_json", cols)
                conn.execute("SELECT 1 FROM publication_versions LIMIT 1")
            finally:
                conn.close()

    def test_migrate_from_v2(self):
        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "library.sqlite"
            # Create a minimal v2-like DB then open_db should migrate
            import sqlite3

            conn = sqlite3.connect(str(db))
            conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
            conn.execute(
                "INSERT INTO meta(key,value) VALUES ('schema_version','2')"
            )
            conn.execute(
                """
                CREATE TABLE publications (
                  id TEXT PRIMARY KEY, doi TEXT, pmid TEXT, pmcid TEXT, arxiv TEXT,
                  openalex TEXT, title TEXT, abstract TEXT, authors_json TEXT,
                  year INTEGER, venue TEXT, type TEXT, tags_json TEXT,
                  needs_review INTEGER NOT NULL DEFAULT 0,
                  created_at TEXT NOT NULL, updated_at TEXT NOT NULL
                )
                """
            )
            conn.execute(
                "CREATE TABLE aliases (alias TEXT PRIMARY KEY, publication_id TEXT NOT NULL)"
            )
            conn.execute(
                """
                CREATE TABLE collections (
                  id TEXT PRIMARY KEY, name TEXT NOT NULL, description TEXT,
                  created_at TEXT NOT NULL, path TEXT NOT NULL
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE collection_pubs (
                  collection_id TEXT NOT NULL, publication_id TEXT NOT NULL,
                  added_at TEXT NOT NULL, source_run_id TEXT,
                  PRIMARY KEY (collection_id, publication_id)
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE runs (
                  id TEXT PRIMARY KEY, collection_id TEXT NOT NULL, query TEXT,
                  sources_json TEXT, config_json TEXT, started_at TEXT,
                  finished_at TEXT, hits_path TEXT
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE run_hits (
                  id INTEGER PRIMARY KEY AUTOINCREMENT, run_id TEXT NOT NULL,
                  publication_id TEXT, raw_json TEXT, source TEXT,
                  fetched_at TEXT, dedupe_action TEXT NOT NULL
                )
                """
            )
            conn.execute(
                "CREATE TABLE theme_nodes (id TEXT PRIMARY KEY, collection_id TEXT NOT NULL, parent_id TEXT, name TEXT NOT NULL, sort INTEGER NOT NULL DEFAULT 0)"
            )
            conn.execute(
                "CREATE TABLE theme_pubs (node_id TEXT NOT NULL, publication_id TEXT NOT NULL, PRIMARY KEY (node_id, publication_id))"
            )
            conn.commit()
            conn.close()
            conn2 = open_db(db)
            try:
                v = conn2.execute(
                    "SELECT value FROM meta WHERE key='schema_version'"
                ).fetchone()["value"]
                self.assertEqual(v, str(SCHEMA_VERSION))
                cols = {r[1] for r in conn2.execute("PRAGMA table_info(publications)")}
                self.assertIn("review_reasons_json", cols)
                conn2.execute("SELECT id FROM publication_versions LIMIT 1")
            finally:
                conn2.close()


class HumanQueueTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db = Path(self.tmp.name) / "library.sqlite"
        init_db(self.db)
        self.conn = open_db(self.db)
        now = utc_now()
        self.conn.execute(
            "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
            ("col-a", "A", "", now, "collections/col-a"),
        )
        self.conn.commit()

    def tearDown(self):
        self.conn.close()
        self.tmp.cleanup()

    def test_empty_queue(self):
        out = list_human_queue(self.conn)
        self.assertTrue(out["ok"])
        self.assertEqual(out["n"], 0)
        self.assertEqual(out["items"], [])

    def test_soft_conflict_and_shell_in_queue(self):
        r1 = resolve_or_create(
            self.conn,
            {
                "doi": "10.1000/soft.a",
                "title": "Soft Fingerprint Collision Example About Exhaustion",
                "authors": ["Doe, Jane"],
                "year": 2019,
            },
        )
        r2 = resolve_or_create(
            self.conn,
            {
                "doi": "10.1000/soft.b",
                "title": "Soft Fingerprint Collision Example About Exhaustion",
                "authors": ["Doe, Jane"],
                "year": 2019,
            },
        )
        shell = resolve_or_create(
            self.conn,
            {"doi": "10.1000/shell.only", "title": "10.1000/shell.only"},
        )
        self.conn.commit()
        self.assertEqual(r2["action"], "needs_review")
        link_to_collection(self.conn, "col-a", r1["publication_id"])
        link_to_collection(self.conn, "col-a", r2["publication_id"])
        link_to_collection(self.conn, "col-a", shell["publication_id"])
        self.conn.commit()

        out = list_human_queue(self.conn, collection_id="col-a")
        self.assertGreaterEqual(out["n"], 2)
        reasons = {i["publication_id"]: i["reasons"] for i in out["items"]}
        self.assertIn("soft_fingerprint_conflict", reasons[r2["publication_id"]])
        self.assertIn("empty_unenriched_shell", reasons[shell["publication_id"]])

        # Filter: unknown collection → empty
        empty = list_human_queue(self.conn, collection_id="col-none")
        self.assertEqual(empty["n"], 0)

    def test_conflict_reason_stored(self):
        resolve_or_create(self.conn, {"doi": "10.1000/a", "title": "A", "year": 2021, "authors": ["X"]})
        resolve_or_create(self.conn, {"pmid": "222", "title": "B", "year": 2022, "authors": ["Y"]})
        self.conn.commit()
        r = resolve_or_create(
            self.conn,
            {"doi": "10.1000/a", "pmid": "222", "title": "Conflict", "year": 2023, "authors": ["Z"]},
        )
        self.conn.commit()
        self.assertEqual(r["action"], "needs_review")
        row = self.conn.execute(
            "SELECT review_reasons_json FROM publications WHERE id=?",
            (r["publication_id"],),
        ).fetchone()
        self.assertIn("conflicting_strong_ids", row["review_reasons_json"])


class CoverageTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db = Path(self.tmp.name) / "library.sqlite"
        init_db(self.db)
        self.conn = open_db(self.db)
        now = utc_now()
        self.conn.execute(
            "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
            ("col-c", "C", "", now, "collections/col-c"),
        )
        self.conn.execute(
            "INSERT INTO runs(id,collection_id,query,sources_json,started_at,finished_at) VALUES (?,?,?,?,?,?)",
            ("run-1", "col-c", "exhaustion", '["pubmed"]', now, now),
        )
        self.conn.commit()

    def tearDown(self):
        self.conn.close()
        self.tmp.cleanup()

    def test_coverage_counts(self):
        a = resolve_or_create(
            self.conn,
            {
                "doi": "10.1000/cov.1",
                "title": "Coverage Paper One",
                "abstract": "Has abstract",
                "year": 2020,
                "authors": ["A"],
            },
        )
        b = resolve_or_create(
            self.conn,
            {"doi": "10.1000/cov.shell", "title": "10.1000/cov.shell"},
        )
        self.conn.commit()
        link_to_collection(self.conn, "col-c", a["publication_id"], source_run_id="run-1")
        link_to_collection(self.conn, "col-c", b["publication_id"], source_run_id="run-1")
        record_run_hit(
            self.conn,
            run_id="run-1",
            publication_id=a["publication_id"],
            raw={"doi": "10.1000/cov.1"},
            source="pubmed",
            dedupe_action="inserted",
        )
        record_run_hit(
            self.conn,
            run_id="run-1",
            publication_id=a["publication_id"],
            raw={"doi": "10.1000/cov.1"},
            source="europepmc",
            dedupe_action="merged",
        )
        record_run_hit(
            self.conn,
            run_id="run-1",
            publication_id=b["publication_id"],
            raw={"doi": "10.1000/cov.shell"},
            source="pubmed",
            dedupe_action="inserted",
        )
        self.conn.commit()

        rep = coverage_report(self.conn, run_id="run-1")
        self.assertTrue(rep["ok"])
        self.assertEqual(rep["hits"]["n"], 3)
        self.assertEqual(rep["hits"]["unique_publication_ids"], 2)
        self.assertEqual(rep["hits"]["dedupe_action_counts"]["inserted"], 2)
        self.assertEqual(rep["hits"]["dedupe_action_counts"]["merged"], 1)
        self.assertEqual(rep["collection"]["n_publications"], 2)
        self.assertGreaterEqual(rep["collection"]["n_unenriched_shells"], 1)
        self.assertGreaterEqual(rep["collection"]["n_missing_abstract"], 1)


class VersionChainTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db = Path(self.tmp.name) / "library.sqlite"
        init_db(self.db)
        self.conn = open_db(self.db)

    def tearDown(self):
        self.conn.close()
        self.tmp.cleanup()

    def _seed_pair(self):
        pre = resolve_or_create(
            self.conn,
            {
                "doi": "10.1101/2020.01.01.123456",
                "arxiv": "2020.01456",
                "title": "Amazing Exhaustion Study In Mice And Humans",
                "authors": [{"family": "Smith"}],
                "year": 2020,
                "abstract": "Preprint abstract about exhaustion.",
                "type": "preprint",
            },
        )
        jou = resolve_or_create(
            self.conn,
            {
                "doi": "10.1038/s41586-020-0001-x",
                "pmid": "333111",
                "title": "Amazing Exhaustion Study In Mice And Humans",
                "authors": [{"family": "Smith"}],
                "year": 2020,
                "abstract": "Journal version with more text about exhaustion and LAYN.",
                "type": "journal-article",
            },
        )
        self.conn.commit()
        return pre, jou

    def test_detect_and_propose(self):
        pre, jou = self._seed_pair()
        props = detect_version_candidates(self.conn)
        self.assertGreaterEqual(len(props), 1)
        hit = next(
            p
            for p in props
            if p["from_id"] == pre["publication_id"] and p["to_id"] == jou["publication_id"]
        )
        self.assertTrue(hit["high_confidence"])
        self.assertEqual(hit["relation"], "preprint_of")

        link = propose_version_link(
            self.conn,
            pre["publication_id"],
            jou["publication_id"],
            relation="preprint_of",
            confidence=hit["confidence"],
            needs_confirm=1,
            evidence=hit["evidence"],
        )
        self.conn.commit()
        self.assertTrue(link["ok"])
        q = list_human_queue(self.conn)
        self.assertGreaterEqual(q["version_links_pending"], 1)
        self.assertTrue(
            any("pending_version_link" in i["reasons"] for i in q["items"])
        )

    def test_confirm_merges_to_journal_doi(self):
        pre, jou = self._seed_pair()
        link = propose_version_link(
            self.conn,
            pre["publication_id"],
            jou["publication_id"],
            relation="preprint_of",
            confidence=0.9,
            needs_confirm=1,
            evidence={"test": True},
        )
        self.conn.commit()
        out = confirm_version_link(self.conn, link["link_id"], merge=True)
        self.conn.commit()
        self.assertTrue(out["ok"])
        self.assertTrue(out["merged"])
        keeper = out["merge"]["publication_id"]
        self.assertEqual(keeper, "doi:10.1038/s41586-020-0001-x")
        n = self.conn.execute("SELECT COUNT(*) c FROM publications").fetchone()["c"]
        self.assertEqual(n, 1)
        row = self.conn.execute(
            "SELECT doi, pmid, arxiv, needs_review FROM publications WHERE id=?",
            (keeper,),
        ).fetchone()
        self.assertEqual(row["doi"], "10.1038/s41586-020-0001-x")
        self.assertEqual(row["pmid"], "333111")
        # preprint DOI kept as alias
        alias = self.conn.execute(
            "SELECT publication_id FROM aliases WHERE alias=?",
            ("doi:10.1101/2020.01.01.123456",),
        ).fetchone()
        self.assertEqual(alias["publication_id"], keeper)
        ax = self.conn.execute(
            "SELECT publication_id FROM aliases WHERE alias=?",
            ("arxiv:2020.01456",),
        ).fetchone()
        self.assertEqual(ax["publication_id"], keeper)

    def test_apply_auto_merge(self):
        self._seed_pair()
        out = apply_version_detection(self.conn, auto_merge_high_confidence=True)
        self.conn.commit()
        self.assertTrue(out["ok"])
        self.assertGreaterEqual(out["n_merged"], 1)
        n = self.conn.execute("SELECT COUNT(*) c FROM publications").fetchone()["c"]
        self.assertEqual(n, 1)

    def test_merge_preprint_journal_direct(self):
        pre, jou = self._seed_pair()
        m = merge_preprint_journal(
            self.conn, pre["publication_id"], jou["publication_id"], evidence={"direct": True}
        )
        self.conn.commit()
        self.assertTrue(m["ok"])
        self.assertEqual(m["publication_id"], "doi:10.1038/s41586-020-0001-x")


class CliSmokeTests(unittest.TestCase):
    def test_cli_scripts(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            conn = open_db(db)
            now = utc_now()
            conn.execute(
                "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
                ("col-cli", "CLI", "", now, "collections/col-cli"),
            )
            conn.execute(
                "INSERT INTO runs(id,collection_id,query,started_at) VALUES (?,?,?,?)",
                ("run-cli", "col-cli", "q", now),
            )
            r = resolve_or_create(
                conn,
                {"doi": "10.1000/cli.1", "title": "CLI Paper", "year": 2021, "authors": ["Z"], "abstract": "x"},
            )
            link_to_collection(conn, "col-cli", r["publication_id"], source_run_id="run-cli")
            record_run_hit(
                conn,
                run_id="run-cli",
                publication_id=r["publication_id"],
                raw={},
                source="pubmed",
                dedupe_action="inserted",
            )
            conn.commit()
            conn.close()

            py = sys.executable
            rq = subprocess.check_output(
                [
                    py,
                    str(ROOT / "skills/literature-library/scripts/review_queue.py"),
                    "--project_root",
                    str(root),
                    "--collection_id",
                    "col-cli",
                ],
                text=True,
            )
            self.assertTrue(json.loads(rq)["ok"])

            cov = subprocess.check_output(
                [
                    py,
                    str(ROOT / "skills/literature-library/scripts/coverage.py"),
                    "--project_root",
                    str(root),
                    "--run_id",
                    "run-cli",
                ],
                text=True,
            )
            cov_j = json.loads(cov)
            self.assertTrue(cov_j["ok"])
            self.assertEqual(cov_j["hits"]["n"], 1)

            ver = subprocess.check_output(
                [
                    py,
                    str(ROOT / "skills/literature-library/scripts/versions.py"),
                    "--project_root",
                    str(root),
                    "--action",
                    "detect",
                ],
                text=True,
            )
            self.assertTrue(json.loads(ver)["ok"])


if __name__ == "__main__":
    unittest.main()
