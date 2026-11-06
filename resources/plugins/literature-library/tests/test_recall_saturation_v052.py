#!/usr/bin/env python3
"""Offline tests for 0.5.2: gold-set recall, OOW leakage, uncontacted sources, snowball saturation."""
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
    PROTOCOL_VERSION,
    assess_snowball_saturation,
    coverage_report,
    evaluate_gold_set_recall,
    evaluate_out_of_window_controls,
    evaluate_uncontacted_sources,
    init_db,
    link_to_collection,
    open_db,
    record_run_hit,
    resolve_or_create,
    utc_now,
)


def _run(cmd: list[str]) -> dict:
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode not in (0, 2):
        raise AssertionError(f"rc={r.returncode}\nstdout={r.stdout}\nstderr={r.stderr}")
    line = (r.stdout or "").strip().splitlines()[-1]
    return json.loads(line), r.returncode


class ProtocolGoldFieldsTests(unittest.TestCase):
    def test_run_config_freezes_gold_and_oow(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            conn = open_db(db)
            col_path = root / "collections" / "col-g"
            col_path.mkdir(parents=True)
            conn.execute(
                "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
                ("col-g", "G", "", utc_now(), str(col_path)),
            )
            conn.commit()
            conn.close()
            out, rc = _run(
                [
                    sys.executable,
                    str(ROOT / "skills/literature-library/scripts/run_config.py"),
                    "--project_root",
                    str(root),
                    "--collection_id",
                    "col-g",
                    "--query",
                    "exhaustion",
                    "--sources",
                    "pubmed,openalex,doaj",
                    "--date_from",
                    "2018-01-01",
                    "--date_to",
                    "2024-12-31",
                    "--gold_set",
                    "10.1000/gold.1,pmid:1001,10.1000/gold.2,10.1000/gold.3",
                    "--out_of_window_controls",
                    "10.1000/oow.old",
                    "--run_id",
                    "run-g1",
                ]
            )
            self.assertEqual(rc, 0)
            cfg = out["config"]
            self.assertEqual(cfg["protocol_version"], "0.5.2")
            self.assertEqual(PROTOCOL_VERSION, "0.5.2")
            self.assertEqual(len(cfg["gold_set"]), 4)
            self.assertTrue(cfg["gold_set"][0].startswith("doi:"))
            self.assertEqual(cfg["out_of_window_controls"], ["doi:10.1000/oow.old"])
            disk = json.loads((col_path / "runs/run-g1/run.json").read_text())
            self.assertEqual(disk["gold_set"], cfg["gold_set"])


class GoldRecallTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db = Path(self.tmp.name) / "library.sqlite"
        init_db(self.db)
        self.conn = open_db(self.db)
        now = utc_now()
        self.col_path = Path(self.tmp.name) / "collections" / "col-r"
        self.col_path.mkdir(parents=True)
        self.conn.execute(
            "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
            ("col-r", "R", "", now, str(self.col_path)),
        )
        self.conn.commit()

    def tearDown(self):
        self.conn.close()
        self.tmp.cleanup()

    def test_found_missed_and_evidence_weak(self):
        a = resolve_or_create(
            self.conn,
            {"doi": "10.1000/gold.1", "title": "Gold One Found", "year": 2019, "authors": ["A"]},
        )
        resolve_or_create(
            self.conn,
            {"doi": "10.1000/other", "title": "Other", "year": 2020, "authors": ["B"]},
        )
        link_to_collection(self.conn, "col-r", a["publication_id"], source_run_id="run-r")
        run_dir = self.col_path / "runs" / "run-r"
        run_dir.mkdir(parents=True)
        cfg = {
            "protocol_version": "0.5.2",
            "run_id": "run-r",
            "collection_id": "col-r",
            "query": "q",
            "sources": ["pubmed", "openalex", "doaj"],
            "sources_executed": ["pubmed"],
            "providers": [
                {"source": "pubmed", "provider": "pubmed", "n": 1},
                {"source": "openalex", "provider": "error", "error": "timeout"},
            ],
            "gold_set": [
                "doi:10.1000/gold.1",
                "doi:10.1000/gold.missing",
                "pmid:999001",
            ],
            "out_of_window_controls": ["doi:10.1000/oow.control"],
            "date_from": "2018-01-01",
            "limit": 10,
        }
        (run_dir / "run.json").write_text(json.dumps(cfg), encoding="utf-8")
        self.conn.execute(
            """
            INSERT INTO runs(id,collection_id,query,sources_json,config_json,started_at,hits_path,status)
            VALUES (?,?,?,?,?,?,?,?)
            """,
            (
                "run-r",
                "col-r",
                "q",
                json.dumps(["pubmed", "openalex", "doaj"]),
                json.dumps(cfg),
                utc_now(),
                str(run_dir / "hits.jsonl"),
                "searched",
            ),
        )
        record_run_hit(
            self.conn,
            run_id="run-r",
            publication_id=a["publication_id"],
            raw={"doi": "10.1000/gold.1"},
            source="pubmed",
            dedupe_action="inserted",
        )
        self.conn.commit()

        rep = coverage_report(self.conn, run_id="run-r")
        self.assertTrue(rep["ok"])
        gold = rep["recall_evidence"]["gold_set_recall"]
        self.assertEqual(gold["gold_set_n"], 3)
        self.assertTrue(gold["evidence_weak"])  # <5
        self.assertTrue(rep["evidence_weak"])
        self.assertEqual(gold["n_found"], 1)
        self.assertEqual(gold["n_missed"], 2)
        found_ids = {x["id"] for x in gold["found"]}
        missed_ids = {x["id"] for x in gold["missed"]}
        self.assertIn("doi:10.1000/gold.1", found_ids)
        self.assertIn("doi:10.1000/gold.missing", missed_ids)

        # empty gold → explicit
        empty = evaluate_gold_set_recall(self.conn, [], corpus_pub_ids={a["publication_id"]})
        self.assertTrue(empty["gold_set_empty"])
        self.assertTrue(empty["evidence_weak"])
        self.assertIn("empty", empty["note"].lower())

    def test_oow_leakage_detected(self):
        leaked = resolve_or_create(
            self.conn,
            {
                "doi": "10.1000/oow.leak",
                "title": "Should Be Outside Window",
                "year": 2010,
                "authors": ["Z"],
            },
        )
        link_to_collection(self.conn, "col-r", leaked["publication_id"], source_run_id="run-oow")
        run_dir = self.col_path / "runs" / "run-oow"
        run_dir.mkdir(parents=True)
        cfg = {
            "sources": ["pubmed"],
            "sources_executed": ["pubmed"],
            "providers": [{"source": "pubmed", "provider": "pubmed"}],
            "gold_set": [],
            "out_of_window_controls": ["doi:10.1000/oow.leak", "doi:10.1000/oow.absent"],
        }
        (run_dir / "run.json").write_text(json.dumps(cfg), encoding="utf-8")
        self.conn.execute(
            """
            INSERT INTO runs(id,collection_id,query,sources_json,config_json,started_at,hits_path,status)
            VALUES (?,?,?,?,?,?,?,?)
            """,
            (
                "run-oow",
                "col-r",
                "q",
                '["pubmed"]',
                json.dumps(cfg),
                utc_now(),
                str(run_dir / "hits.jsonl"),
                "searched",
            ),
        )
        record_run_hit(
            self.conn,
            run_id="run-oow",
            publication_id=leaked["publication_id"],
            raw={"doi": "10.1000/oow.leak"},
            source="pubmed",
            dedupe_action="inserted",
        )
        self.conn.commit()
        rep = coverage_report(self.conn, run_id="run-oow")
        oow = rep["recall_evidence"]["out_of_window"]
        self.assertEqual(oow["n_leaked"], 1)
        self.assertFalse(oow["filter_ok"])
        self.assertEqual(oow["leaked"][0]["id"], "doi:10.1000/oow.leak")
        self.assertEqual(oow["n_absent_ok"], 1)

        ok = evaluate_out_of_window_controls(
            self.conn, ["doi:10.1000/oow.absent"], corpus_pub_ids=set()
        )
        self.assertTrue(ok["filter_ok"])

    def test_uncontacted_sources_listed(self):
        out = evaluate_uncontacted_sources(
            configured_sources=["pubmed", "openalex", "doaj"],
            sources_executed=["pubmed"],
            providers=[
                {"source": "pubmed", "provider": "pubmed", "n": 3},
                {"source": "openalex", "provider": "error", "error": "boom"},
            ],
            hit_source_counts={"pubmed": 3},
        )
        reasons = {i["source"]: i["reason"] for i in out["uncontacted"]}
        self.assertIn("doaj", reasons)
        self.assertEqual(reasons["doaj"], "not_in_sources_executed")
        self.assertIn("openalex", reasons)
        # section never blank-implies-ok
        empty_cfg = evaluate_uncontacted_sources(configured_sources=[])
        self.assertIn("no sources", empty_cfg["note"].lower())


class SaturationTests(unittest.TestCase):
    def test_assess_stop_reasons(self):
        sat = assess_snowball_saturation(
            [
                {"round": 1, "n_seeds": 10, "n_new_included": 8},
                {"round": 2, "n_seeds": 8, "n_new_included": 6},
            ],
            max_rounds=2,
        )
        self.assertFalse(sat["saturated"])
        self.assertEqual(sat["stop_reason"], "max_rounds_not_saturated")
        self.assertFalse(sat["closure_claimed"])

        sat2 = assess_snowball_saturation(
            [{"round": 1, "n_seeds": 5, "n_new_included": 0}],
            max_rounds=5,
        )
        self.assertTrue(sat2["saturated"])
        self.assertEqual(sat2["stop_reason"], "no_new_eligible_seeds")

        sat3 = assess_snowball_saturation(
            [
                {"round": 1, "n_seeds": 100, "n_new_included": 2},
                {"round": 2, "n_seeds": 100, "n_new_included": 3},
            ],
            max_rounds=5,
        )
        self.assertTrue(sat3["saturated"])
        self.assertEqual(sat3["stop_reason"], "consecutive_low_yield")

    def test_expand_cli_mock_saturation_fields(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            conn = open_db(db)
            col_path = root / "collections" / "col-e"
            col_path.mkdir(parents=True)
            conn.execute(
                "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
                ("col-e", "E", "", utc_now(), str(col_path)),
            )
            conn.commit()
            conn.close()
            # max_rounds=1 with mock hits → NOT saturated (closure not claimed)
            out1, rc1 = _run(
                [
                    sys.executable,
                    str(ROOT / "skills/literature-library/scripts/expand.py"),
                    "--project_root",
                    str(root),
                    "--collection_id",
                    "col-e",
                    "--seed_id",
                    "doi:10.1000/seed",
                    "--mock",
                    "--max_rounds",
                    "1",
                ]
            )
            self.assertEqual(rc1, 0)
            self.assertFalse(out1["saturated"])
            self.assertEqual(out1["stop_reason"], "max_rounds_not_saturated")
            self.assertIn("NOT claimed", out1["banner"])

            # max_rounds=3: round1 new hits, round2 same mock → no new → saturated
            out2, rc2 = _run(
                [
                    sys.executable,
                    str(ROOT / "skills/literature-library/scripts/expand.py"),
                    "--project_root",
                    str(root),
                    "--collection_id",
                    "col-e",
                    "--seed_id",
                    "doi:10.1000/seed",
                    "--mock",
                    "--max_rounds",
                    "3",
                ]
            )
            self.assertEqual(rc2, 0)
            self.assertTrue(out2["saturated"])
            self.assertEqual(out2["stop_reason"], "no_new_eligible_seeds")
            self.assertGreaterEqual(len(out2["rounds"]), 2)


if __name__ == "__main__":
    unittest.main()
