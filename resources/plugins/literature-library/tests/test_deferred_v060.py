#!/usr/bin/env python3
"""Offline tests for 0.6.0 deferred items.

schema v6 and PDF text extraction.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "tests" / "fixtures"


from db import SCHEMA_VERSION, init_db, open_db, record_pub_file, utc_now  # noqa: E402
from pdf_extract import extract_pdf_text, extract_with_pdftotext  # noqa: E402
from rebind_files import file_content_hash, rebind_pub_files  # noqa: E402
from ingest import load_records  # noqa: E402


def _col(root: Path, cid: str = "col-v60") -> str:
    db = root / "library.sqlite"
    init_db(db)
    col_path = root / "collections" / cid
    col_path.mkdir(parents=True)
    conn = open_db(db)
    conn.execute(
        "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
        (cid, "V60", "", utc_now(), str(col_path)),
    )
    conn.commit()
    conn.close()
    return cid


class SchemaV6Tests(unittest.TestCase):
    def test_schema_v6_access_columns(self):
        self.assertEqual(SCHEMA_VERSION, 7)
        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "library.sqlite"
            info = init_db(db)
            self.assertEqual(info["schema_version"], 7)
            conn = open_db(db)
            try:
                cols = {r[1] for r in conn.execute("PRAGMA table_info(pub_files)")}
                self.assertIn("access_origin", cols)
                self.assertIn("access_method", cols)
            finally:
                conn.close()

    def test_migrate_from_v5(self):
        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "library.sqlite"
            # create minimal v5-ish DB then open_db migrates
            import sqlite3

            conn = sqlite3.connect(str(db))
            conn.execute("CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
            conn.execute(
                "INSERT INTO meta(key,value) VALUES ('schema_version','5')"
            )
            conn.execute(
                """
                CREATE TABLE publications (
                  id TEXT PRIMARY KEY, title TEXT, abstract TEXT,
                  authors_json TEXT, year TEXT, venue TEXT, doi TEXT, pmid TEXT,
                  pmcid TEXT, arxiv TEXT, openalex TEXT, url TEXT, type TEXT,
                  needs_review INTEGER DEFAULT 0, created_at TEXT, updated_at TEXT,
                  oa_status TEXT, oa_json TEXT, review_reasons_json TEXT
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE pub_files (
                  id INTEGER PRIMARY KEY AUTOINCREMENT,
                  publication_id TEXT NOT NULL,
                  file_path TEXT NOT NULL,
                  content_hash TEXT,
                  ingest_origin TEXT NOT NULL DEFAULT 'local_file',
                  created_at TEXT NOT NULL,
                  UNIQUE(publication_id, file_path)
                )
                """
            )
            conn.commit()
            conn.close()
            conn2 = open_db(db)
            try:
                row = conn2.execute(
                    "SELECT value FROM meta WHERE key='schema_version'"
                ).fetchone()
                self.assertEqual(row["value"], "7")
                cols = {r[1] for r in conn2.execute("PRAGMA table_info(pub_files)")}
                self.assertIn("access_origin", cols)
                pub_cols = {r[1] for r in conn2.execute("PRAGMA table_info(publications)")}
                self.assertIn("enriched_at", pub_cols)
            finally:
                conn2.close()




class PdfExtractTests(unittest.TestCase):
    def test_pdftotext_on_tiny_fixture(self):
        tiny = FX / "tiny.pdf"
        self.assertTrue(tiny.is_file())
        res = extract_with_pdftotext(tiny)
        self.assertTrue(res["ok"])
        self.assertIn("10.1000/tiny.pdf.001", res["text"])
        self.assertIn("7770001", res["text"])

    def test_load_records_pdf_extracts_ids(self):
        records = load_records(FX / "tiny.pdf", fmt="pdf", try_ocr=False)
        self.assertEqual(len(records), 1)
        r = records[0]
        self.assertEqual(r.get("doi"), "10.1000/tiny.pdf.001")
        self.assertEqual(r.get("pmid"), "7770001")
        self.assertEqual(r.get("extract_method"), "pdftotext")
        # The PDF body is not an abstract; tiny.pdf has no Abstract heading.
        self.assertFalse((r.get("abstract") or "").strip())

    def test_ocr_soft_fail_when_missing(self):
        # Force empty digital text via a shim that prints nothing
        with tempfile.TemporaryDirectory() as td:
            shim_dir = Path(td) / "bin"
            shim_dir.mkdir()
            shim = shim_dir / "pdftotext"
            shim.write_text(
                "#!/bin/sh\n# emit empty\nexit 0\n",
                encoding="utf-8",
            )
            shim.chmod(0o755)
            empty_pdf = Path(td) / "empty.pdf"
            # reuse tiny bytes but empty tool output
            empty_pdf.write_bytes((FX / "tiny.pdf").read_bytes())
            res = extract_pdf_text(
                empty_pdf,
                try_ocr=True,
                pdftotext_bin=str(shim),
                tesseract_bin="/nonexistent/tesseract",
            )
            joined = " ".join(res.get("warnings") or [])
            self.assertTrue(
                "tesseract_missing" in joined or "tesseract_failed" in joined,
                joined,
            )
            # ok may be False when no text — soft, not raise
            self.assertIsInstance(res.get("warnings"), list)








if __name__ == "__main__":
    unittest.main()
