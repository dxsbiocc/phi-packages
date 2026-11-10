#!/usr/bin/env python3
"""Offline tests: runs schema columns and preview_fetch of the opposite preprint/journal record."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "tests" / "fixtures" / "providers"

from db import (  # noqa: E402
    SCHEMA_VERSION,
    init_db,
    link_to_collection,
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
                    str(ROOT / "skills/manage-library/scripts/lib/versions.py"),
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
