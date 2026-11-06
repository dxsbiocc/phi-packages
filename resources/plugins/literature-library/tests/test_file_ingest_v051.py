#!/usr/bin/env python3
"""Offline tests for 0.5.1: file ingest polish + OA verify (no network)."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "tests" / "fixtures"
sys_path_db = ROOT / "skills" / "literature-library" / "scripts"
sys_path_ingest = ROOT / "skills" / "literature-library" / "scripts"

import sys

sys.path.insert(0, str(sys_path_db))
sys.path.insert(0, str(sys_path_ingest))

from db import (  # noqa: E402
    SCHEMA_VERSION,
    init_db,
    open_db,
    utc_now,
)
from oa_verify import (  # noqa: E402
    candidate_oa_urls,
    verify_oa_url,
    verify_record_oa,
)
from ingest import (  # noqa: E402
    extract_ids_from_filename,
    extract_ids_from_text,
    file_content_hash,
    ingest_records,
    load_records,
)


def _col(root: Path, cid: str = "col-f") -> str:
    db = root / "library.sqlite"
    init_db(db)
    col_path = root / "collections" / cid
    col_path.mkdir(parents=True)
    conn = open_db(db)
    conn.execute(
        "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
        (cid, "F", "", utc_now(), str(col_path)),
    )
    conn.commit()
    conn.close()
    return cid


class SchemaV5Tests(unittest.TestCase):
    def test_schema_v5_pub_files_and_oa_cols(self):
        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "library.sqlite"
            info = init_db(db)
            self.assertGreaterEqual(info["schema_version"], 5)
            self.assertGreaterEqual(SCHEMA_VERSION, 5)
            conn = open_db(db)
            try:
                tables = {
                    r[0]
                    for r in conn.execute(
                        "SELECT name FROM sqlite_master WHERE type='table'"
                    )
                }
                self.assertIn("pub_files", tables)
                cols = {r[1] for r in conn.execute("PRAGMA table_info(publications)")}
                self.assertIn("oa_status", cols)
                self.assertIn("oa_json", cols)
            finally:
                conn.close()


class IdExtractionTests(unittest.TestCase):
    def test_extract_from_text_and_filename(self):
        text_ids = extract_ids_from_text(
            "See DOI: 10.1016/j.cell.2017.05.035 and PMID: 28622514 PMC5482351"
        )
        self.assertEqual(text_ids.get("doi"), "10.1016/j.cell.2017.05.035")
        self.assertEqual(text_ids.get("pmid"), "28622514")
        self.assertTrue(str(text_ids.get("pmcid", "")).upper().startswith("PMC"))

        fn = extract_ids_from_filename("10.1016_j.cell.2017.05.035.pdf")
        self.assertEqual(fn.get("doi"), "10.1016/j.cell.2017.05.035")
        fn2 = extract_ids_from_filename("PMID28622514_zheng.pdf")
        self.assertEqual(fn2.get("pmid"), "28622514")


class PdfSidecarIngestTests(unittest.TestCase):
    def test_sidecar_txt_records_file_and_ids(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            src = FX / "pdf_sidecar_zheng.txt"
            records = load_records(src, fmt="pdf_text")
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].get("doi"), "10.1016/j.cell.2017.05.035")
            self.assertEqual(records[0].get("pmid"), "28622514")
            self.assertEqual(records[0].get("ingest_origin"), "local_file")
            out = ingest_records(
                db_path=root / "library.sqlite",
                collection_id=cid,
                records=records,
            )
            self.assertTrue(out["ok"])
            self.assertGreaterEqual(out["n_files_recorded"], 1)
            conn = open_db(root / "library.sqlite")
            try:
                row = conn.execute(
                    "SELECT id, doi, pmid, title FROM publications"
                ).fetchone()
                self.assertIsNotNone(row)
                self.assertEqual(row["doi"], "10.1016/j.cell.2017.05.035")
                files = conn.execute("SELECT * FROM pub_files").fetchall()
                self.assertGreaterEqual(len(files), 1)
                self.assertEqual(files[0]["ingest_origin"], "local_file")
                self.assertTrue(files[0]["content_hash"])
                self.assertEqual(len(files[0]["content_hash"]), 64)
            finally:
                conn.close()

    def test_json_sidecar_filename_doi_and_hash(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            src = FX / "10.1016_j.cell.2017.05.035.json"
            records = load_records(src, fmt="pdf_text")
            self.assertEqual(records[0].get("doi"), "10.1016/j.cell.2017.05.035")
            self.assertEqual(records[0].get("pmid"), "28622514")
            h = file_content_hash(src)
            out = ingest_records(
                db_path=root / "library.sqlite",
                collection_id=cid,
                records=records,
            )
            self.assertTrue(out["ok"])
            self.assertEqual(out["files_recorded"][0]["content_hash"], h)


class RisLocalIngestTests(unittest.TestCase):
    def test_ris_url_doi_extraction(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            src = FX / "sample_local.ris"
            records = load_records(src, fmt="ris")
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].get("doi"), "10.1000/ris.local.002")
            out = ingest_records(
                db_path=root / "library.sqlite",
                collection_id=cid,
                records=records,
            )
            self.assertTrue(out["ok"])
            self.assertGreaterEqual(out["counts"]["inserted"] + out["counts"].get("merged", 0), 1)
            self.assertGreaterEqual(out["n_files_recorded"], 1)


class PdfBinaryStubTests(unittest.TestCase):
    def test_pdf_uses_sidecar_and_registers_pdf_path(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            pdf = root / "paper.pdf"
            pdf.write_bytes(b"%PDF-1.4 fake")
            sidecar = root / "paper.txt"
            sidecar.write_text(
                "---\ntitle: Fake PDF Sidecar\ndoi: 10.1000/fake.pdf.001\n---\nBody PMID: 12345678\n",
                encoding="utf-8",
            )
            records = load_records(pdf, fmt="pdf")
            self.assertEqual(records[0]["doi"], "10.1000/fake.pdf.001")
            self.assertEqual(records[0]["pmid"], "12345678")
            self.assertEqual(Path(records[0]["_local_file"]), pdf.resolve())
            out = ingest_records(
                db_path=root / "library.sqlite",
                collection_id=cid,
                records=records,
            )
            self.assertTrue(out["ok"])
            conn = open_db(root / "library.sqlite")
            try:
                paths = {
                    r["file_path"]
                    for r in conn.execute("SELECT file_path FROM pub_files")
                }
                self.assertTrue(any(str(pdf.resolve()) == p for p in paths))
            finally:
                conn.close()

    def test_pdf_filename_only_shell(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            pdf = root / "10.1000_filename.only.001.pdf"
            pdf.write_bytes(b"%PDF-1.4")
            records = load_records(pdf, fmt="pdf")
            self.assertEqual(records[0].get("doi"), "10.1000/filename.only.001")
            out = ingest_records(
                db_path=root / "library.sqlite",
                collection_id=cid,
                records=records,
            )
            self.assertTrue(out["ok"])
            self.assertEqual(out["n_files_recorded"], 1)


class AutoEnrichLocalShellTests(unittest.TestCase):
    def test_local_shell_auto_enrich_mock(self):
        import subprocess

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            # ID-only sidecar → should auto-enrich with mock
            ids = root / "shell.ids"
            ids.write_text("doi:10.1016/j.cell.2017.05.035\n", encoding="utf-8")
            cmd = [
                sys.executable,
                str(ROOT / "skills/literature-library/scripts/ingest.py"),
                "--project_root",
                str(root),
                "--collection_id",
                cid,
                "--input_path",
                str(ids),
                "--format",
                "ids",
                "--enrich_mock",
            ]
            r = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            out = json.loads(r.stdout.strip().splitlines()[-1])
            self.assertTrue(out.get("ok"))
            self.assertTrue(out.get("enrich_ok"))
            self.assertGreaterEqual((out.get("enrich") or {}).get("n_enriched") or 0, 1)


class OaVerifyTests(unittest.TestCase):
    def test_classify_mock_statuses(self):
        mock_map = json.loads((FX / "oa_mock_head.json").read_text(encoding="utf-8"))
        pdf_url = (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5482351/pdf/nihms862.pdf"
        )
        r1 = verify_oa_url(pdf_url, mock=mock_map[pdf_url])
        self.assertEqual(r1["oa_status"], "verified_oa")
        self.assertFalse(r1.get("stored_pdf_bytes"))
        self.assertTrue(r1.get("mock"))

        r2 = verify_oa_url(
            "https://paywall.example/login?article=1",
            mock=mock_map["https://paywall.example/login?article=1"],
        )
        self.assertEqual(r2["oa_status"], "paywall_or_blocked")

        r3 = verify_oa_url(
            "https://publisher.example/article/123",
            mock=mock_map["https://publisher.example/article/123"],
        )
        self.assertEqual(r3["oa_status"], "not_pdf")

    def test_ingest_verify_oa_stores_json_not_bytes(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            src = FX / "10.1016_j.cell.2017.05.035.json"
            records = load_records(src, fmt="pdf_text")
            mock_map = json.loads((FX / "oa_mock_head.json").read_text(encoding="utf-8"))
            out = ingest_records(
                db_path=root / "library.sqlite",
                collection_id=cid,
                records=records,
                verify_oa=True,
                oa_mock_map=mock_map,
            )
            self.assertTrue(out["ok"])
            self.assertEqual(len(out["oa_checks"]), 1)
            self.assertEqual(out["oa_checks"][0]["oa_status"], "verified_oa")
            self.assertFalse(out["oa_checks"][0].get("stored_pdf_bytes"))
            conn = open_db(root / "library.sqlite")
            try:
                row = conn.execute(
                    "SELECT oa_status, oa_json FROM publications"
                ).fetchone()
                self.assertEqual(row["oa_status"], "verified_oa")
                payload = json.loads(row["oa_json"])
                self.assertEqual(payload["oa_status"], "verified_oa")
                self.assertNotIn("pdf_bytes", payload)
                self.assertFalse(payload.get("stored_pdf_bytes"))
            finally:
                conn.close()

    def test_candidate_urls_from_record(self):
        urls = candidate_oa_urls(
            {
                "url": "https://example.com/a",
                "best_oa_location": {"url": "https://example.com/pdf"},
            }
        )
        self.assertEqual(urls[0], "https://example.com/a")
        self.assertIn("https://example.com/pdf", urls)


class SampleBibStillWorks(unittest.TestCase):
    def test_sample_bib_records_local_file(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            records = load_records(FX / "sample.bib", fmt="bib")
            out = ingest_records(
                db_path=root / "library.sqlite",
                collection_id=cid,
                records=records,
            )
            self.assertTrue(out["ok"])
            self.assertGreaterEqual(out["n_files_recorded"], 1)


if __name__ == "__main__":
    unittest.main()
