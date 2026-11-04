#!/usr/bin/env python3
"""Offline tests for 0.5.0: protocol freeze, coverage gate, dig queue, preview_fetch."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "tests" / "fixtures" / "providers"
sys.path.insert(0, str(ROOT / "skills" / "lit-db" / "scripts"))

from db import (  # noqa: E402
    SCHEMA_VERSION,
    close_run,
    evaluate_coverage_gate,
    init_db,
    link_to_collection,
    list_human_queue,
    open_db,
    preview_fetch_opposite,
    resolve_or_create,
    unresolved_doi_ties,
    utc_now,
)


def _run(cmd: list[str]) -> dict:
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode not in (0, 2):
        raise AssertionError(f"rc={r.returncode}\nstdout={r.stdout}\nstderr={r.stderr}")
    line = (r.stdout or "").strip().splitlines()[-1]
    return json.loads(line), r.returncode


class SchemaV4Tests(unittest.TestCase):
    def test_schema_version_4_runs_gate_cols(self):
        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "library.sqlite"
            info = init_db(db)
            self.assertGreaterEqual(info["schema_version"], 4)
            self.assertGreaterEqual(SCHEMA_VERSION, 4)
            conn = open_db(db)
            try:
                cols = {r[1] for r in conn.execute("PRAGMA table_info(runs)")}
                self.assertIn("status", cols)
                self.assertIn("coverage_path", cols)
                self.assertIn("gate_json", cols)
            finally:
                conn.close()


class ProtocolFreezeCliTests(unittest.TestCase):
    def test_run_config_freezes_protocol_fields(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            conn = open_db(db)
            col_path = root / "collections" / "col-x"
            col_path.mkdir(parents=True)
            conn.execute(
                "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
                ("col-x", "X", "", utc_now(), str(col_path)),
            )
            conn.commit()
            conn.close()
            out, rc = _run(
                [
                    sys.executable,
                    str(ROOT / "skills/lit-run-config/scripts/run_config.py"),
                    "--project_root",
                    str(root),
                    "--collection_id",
                    "col-x",
                    "--query",
                    "T cell exhaustion",
                    "--sources",
                    "pubmed,openalex",
                    "--limit",
                    "7",
                    "--date_from",
                    "2018-01-01",
                    "--run_id",
                    "run-proto1",
                ]
            )
            self.assertEqual(rc, 0)
            self.assertTrue(out["ok"])
            self.assertTrue(out["protocol_frozen"])
            cfg = out["config"]
            self.assertEqual(cfg["protocol_version"], "0.5.2")
            self.assertEqual(cfg["limit_per_source"], 7)
            self.assertIn("provider_notes", cfg)
            self.assertIn("pubmed", cfg["provider_notes"])
            self.assertTrue(cfg["frozen_at"])
            disk = json.loads((col_path / "runs/run-proto1/run.json").read_text())
            self.assertEqual(disk["sources"], ["pubmed", "openalex"])
            self.assertEqual(disk["date_from"], "2018-01-01")

            # mock search writes hits.jsonl + hits.meta.json
            out2, rc2 = _run(
                [
                    sys.executable,
                    str(ROOT / "skills/lit-search/scripts/search.py"),
                    "--project_root",
                    str(root),
                    "--collection_id",
                    "col-x",
                    "--run_id",
                    "run-proto1",
                    "--mock",
                    "--close_run",
                    "--min_unique_pubs",
                    "1",
                ]
            )
            self.assertEqual(rc2, 0)
            self.assertEqual(out2["n_hits"], 2)
            hits = col_path / "runs/run-proto1/hits.jsonl"
            meta = col_path / "runs/run-proto1/hits.meta.json"
            cov = col_path / "runs/run-proto1/coverage.json"
            self.assertTrue(hits.is_file())
            self.assertTrue(meta.is_file())
            self.assertTrue(cov.is_file())
            self.assertEqual(out2["run_status"], "passed")
            run_row_cfg = json.loads((col_path / "runs/run-proto1/run.json").read_text())
            self.assertEqual(run_row_cfg["status"], "passed")


class CoverageGateTests(unittest.TestCase):
    def test_evaluate_and_close_run_fail(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            conn = open_db(db)
            col_path = root / "collections" / "col-g"
            run_dir = col_path / "runs" / "run-g"
            run_dir.mkdir(parents=True)
            now = utc_now()
            conn.execute(
                "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
                ("col-g", "G", "", now, str(col_path)),
            )
            cfg = {"run_id": "run-g", "query": "q", "sources": ["pubmed"], "limit": 5}
            (run_dir / "run.json").write_text(json.dumps(cfg), encoding="utf-8")
            (run_dir / "hits.jsonl").write_text("", encoding="utf-8")
            conn.execute(
                """
                INSERT INTO runs(id,collection_id,query,sources_json,config_json,started_at,hits_path,status)
                VALUES (?,?,?,?,?,?,?,?)
                """,
                (
                    "run-g",
                    "col-g",
                    "q",
                    '["pubmed"]',
                    json.dumps(cfg),
                    now,
                    str(run_dir / "hits.jsonl"),
                    "searched",
                ),
            )
            conn.commit()
            gate = evaluate_coverage_gate(
                {"ok": True, "hits": {"n": 0, "unique_publication_ids": 0}, "collection": {"n_needs_review": 0, "n_version_links_pending": 0, "n_unenriched_shells": 0}},
                min_unique_pubs=3,
                fail_under_coverage=True,
            )
            self.assertFalse(gate["passed"])
            self.assertTrue(gate["should_fail_cli"])
            out = close_run(
                conn,
                "run-g",
                min_unique_pubs=3,
                fail_under_coverage=True,
            )
            conn.commit()
            self.assertEqual(out["status"], "needs_attention")
            self.assertTrue(out["should_fail_cli"])
            self.assertTrue(Path(out["coverage_json"]).is_file())
            status = conn.execute("SELECT status FROM runs WHERE id='run-g'").fetchone()["status"]
            self.assertEqual(status, "needs_attention")
            conn.close()

            out_cli, rc = _run(
                [
                    sys.executable,
                    str(ROOT / "skills/lit-run-close/scripts/run_close.py"),
                    "--project_root",
                    str(root),
                    "--run_id",
                    "run-g",
                    "--min_unique_pubs",
                    "3",
                    "--fail_under_coverage",
                ]
            )
            self.assertEqual(rc, 2)
            self.assertEqual(out_cli["status"], "needs_attention")


class DigQueueTests(unittest.TestCase):
    def test_dig_silent_when_empty(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            out, rc = _run(
                [
                    sys.executable,
                    str(ROOT / "skills/lit-review-queue/scripts/review_queue.py"),
                    "--project_root",
                    str(root),
                    "--dig",
                ]
            )
            self.assertEqual(rc, 0)
            self.assertTrue(out["silent"])
            self.assertEqual(out["n"], 0)
            self.assertEqual(out["items"], [])

    def test_dig_with_items_and_write_report(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            conn = open_db(db)
            now = utc_now()
            conn.execute(
                "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
                ("col-d", "D", "", now, str(root / "collections/col-d")),
            )
            r = resolve_or_create(
                conn,
                {
                    "doi": "10.1000/shell.only",
                    "title": "[shell] doi:10.1000/shell.only",
                    "abstract": "",
                },
            )
            link_to_collection(conn, "col-d", r["publication_id"])
            # force needs_review soft
            conn.execute(
                "UPDATE publications SET needs_review=1, review_reasons_json=? WHERE id=?",
                (json.dumps(["soft_fingerprint_conflict"]), r["publication_id"]),
            )
            conn.commit()
            q = list_human_queue(conn, collection_id="col-d")
            self.assertGreater(q["n"], 0)
            self.assertTrue(q["items"][0].get("suggested_next_action"))
            conn.close()

            report = root / "dig.md"
            out, rc = _run(
                [
                    sys.executable,
                    str(ROOT / "skills/lit-review-queue/scripts/review_queue.py"),
                    "--project_root",
                    str(root),
                    "--collection_id",
                    "col-d",
                    "--dig",
                    "--write_report",
                    str(report),
                ]
            )
            self.assertEqual(rc, 0)
            self.assertFalse(out.get("silent"))
            self.assertGreater(out["n"], 0)
            self.assertTrue(report.is_file())
            self.assertIn("suggested_next_action", out["items"][0])


class PreviewFetchTests(unittest.TestCase):
    def test_unresolved_and_preview_fetch_opposite(self):
        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "library.sqlite"
            init_db(db)
            conn = open_db(db)
            pre = resolve_or_create(
                conn,
                {
                    "doi": "10.1101/2020.01.01.123456",
                    "title": "Amazing Exhaustion Study In Mice And Humans",
                    "authors": [{"family": "Smith"}],
                    "year": 2020,
                    "type": "preprint",
                },
            )
            conn.commit()
            ties = [
                {
                    "from_doi": "10.1101/2020.01.01.123456",
                    "to_doi": "10.1038/s41586-020-0001-x",
                    "relation": "preprint_of",
                    "confidence": 0.95,
                    "high_confidence": True,
                    "provider": "mock",
                }
            ]
            unresolved = unresolved_doi_ties(conn, ties)
            self.assertEqual(len(unresolved), 1)
            self.assertEqual(unresolved[0]["missing_doi"], "10.1038/s41586-020-0001-x")
            self.assertEqual(unresolved[0]["present_id"], pre["publication_id"])

            hit = json.loads((FX / "preview_fetch_dois.json").read_text())["dois"][
                "10.1038/s41586-020-0001-x"
            ]
            pr = preview_fetch_opposite(
                conn,
                missing_doi="10.1038/s41586-020-0001-x",
                present_id=pre["publication_id"],
                relation="preprint_of",
                confidence=0.95,
                hit=hit,
                missing_role="to",
            )
            conn.commit()
            self.assertTrue(pr["ok"])
            self.assertTrue(pr["metadata_only"])
            new_id = pr["publication_id"]
            row = conn.execute("SELECT * FROM publications WHERE id=?", (new_id,)).fetchone()
            self.assertEqual(row["doi"], "10.1038/s41586-020-0001-x")
            self.assertEqual(int(row["needs_review"]), 1)
            reasons = json.loads(row["review_reasons_json"])
            self.assertIn("preview_fetched_opposite", reasons)
            self.assertIn("pending_version_link", reasons)
            n_links = conn.execute(
                "SELECT COUNT(*) c FROM publication_versions WHERE needs_confirm=1"
            ).fetchone()["c"]
            self.assertGreaterEqual(n_links, 1)
            conn.close()

    def test_versions_cli_preview_fetch_mock(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            conn = open_db(db)
            col_path = root / "collections" / "col-v"
            col_path.mkdir(parents=True)
            now = utc_now()
            conn.execute(
                "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
                ("col-v", "V", "", now, str(col_path)),
            )
            pre = resolve_or_create(
                conn,
                {
                    "doi": "10.1101/2020.01.01.123456",
                    "title": "Amazing Exhaustion Study In Mice And Humans",
                    "authors": [{"family": "Smith"}],
                    "year": 2020,
                    "type": "preprint",
                },
            )
            link_to_collection(conn, "col-v", pre["publication_id"])
            conn.commit()
            conn.close()

            out, rc = _run(
                [
                    sys.executable,
                    str(ROOT / "skills/lit-versions/scripts/versions.py"),
                    "--project_root",
                    str(root),
                    "--collection_id",
                    "col-v",
                    "--action",
                    "apply",
                    "--mock",
                    "--preview_fetch",
                    "--fixture_dir",
                    str(FX),
                ]
            )
            self.assertEqual(rc, 0)
            self.assertTrue(out["ok"])
            self.assertTrue(out.get("preview_fetch"))
            self.assertGreaterEqual(out.get("n_unresolved", 0) + len(out.get("preview_results") or []), 1)
            # opposite journal should now exist
            conn = open_db(db)
            n = conn.execute(
                "SELECT COUNT(*) c FROM publications WHERE doi=?",
                ("10.1038/s41586-020-0001-x",),
            ).fetchone()["c"]
            self.assertEqual(n, 1)
            conn.close()


if __name__ == "__main__":
    unittest.main()
