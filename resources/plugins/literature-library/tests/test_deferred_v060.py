#!/usr/bin/env python3
"""Offline tests for 0.6.0 deferred items.

path-watch / PDF extract / institutional / MCP mock / research-evidence bridge.
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
REEFV = ROOT.parent / "research-evidence"

sys.path.insert(0, str(ROOT / "skills" / "literature-library" / "scripts"))

from db import SCHEMA_VERSION, init_db, open_db, record_pub_file, utc_now  # noqa: E402
from pdf_extract import extract_pdf_text, extract_with_pdftotext  # noqa: E402
from rebind_files import file_content_hash, rebind_pub_files  # noqa: E402
from watch_rebind import watch_rebind  # noqa: E402
from ingest import load_records  # noqa: E402
from providers_mcp import (  # noqa: E402
    related_publications_mcp,
    resolve_backend,
    search_publications_mcp,
)
from institutional_access import (  # noqa: E402
    institutional_fetch_and_mount,
    mount_authorized_file,
)


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
        self.assertEqual(SCHEMA_VERSION, 6)
        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "library.sqlite"
            info = init_db(db)
            self.assertEqual(info["schema_version"], 6)
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
                self.assertEqual(row["value"], "6")
                cols = {r[1] for r in conn2.execute("PRAGMA table_info(pub_files)")}
                self.assertIn("access_origin", cols)
            finally:
                conn2.close()


class WatchRebindTests(unittest.TestCase):
    def test_watch_once_dry_and_confirm(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            old_dir = root / "old"
            new_dir = root / "new"
            old_dir.mkdir()
            new_dir.mkdir()
            payload = b"%PDF-1.4 watch-rebind-bytes-v060"
            old_path = old_dir / "paper.pdf"
            old_path.write_bytes(payload)
            digest = file_content_hash(old_path)

            conn = open_db(root / "library.sqlite")
            conn.execute(
                "INSERT INTO publications(id,title,created_at,updated_at) VALUES (?,?,?,?)",
                ("doi:10.watch/1", "Watch paper", utc_now(), utc_now()),
            )
            conn.execute(
                "INSERT INTO collection_pubs(collection_id,publication_id,added_at) VALUES (?,?,?)",
                (cid, "doi:10.watch/1", utc_now()),
            )
            record_pub_file(
                conn,
                publication_id="doi:10.watch/1",
                file_path=str(old_path),
                content_hash=digest,
                ingest_origin="local_file",
            )
            conn.commit()
            conn.close()

            # move file
            new_path = new_dir / "paper-moved.pdf"
            shutil.move(str(old_path), str(new_path))

            dry = watch_rebind(
                project_root=root,
                scan_dir=new_dir,
                collection_id=cid,
                confirm=False,
                once=True,
            )
            self.assertTrue(dry["ok"])
            self.assertTrue(dry["dry_run"])
            self.assertEqual(dry["n_passes"], 1)
            self.assertGreaterEqual(dry["last"]["n_matched"], 1)

            done = watch_rebind(
                project_root=root,
                scan_dir=new_dir,
                collection_id=cid,
                confirm=True,
                once=True,
            )
            self.assertTrue(done["ok"])
            self.assertFalse(done["dry_run"])
            self.assertGreaterEqual(done["last"]["n_updated"], 1)

            conn = open_db(root / "library.sqlite")
            row = conn.execute(
                "SELECT file_path FROM pub_files WHERE publication_id=?",
                ("doi:10.watch/1",),
            ).fetchone()
            conn.close()
            self.assertEqual(row["file_path"], str(new_path.resolve()))


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
        self.assertTrue((r.get("abstract") or "").strip())

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


class InstitutionalTests(unittest.TestCase):
    def test_mount_user_provided_file(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            pdf = root / "auth.pdf"
            shutil.copy(FX / "tiny.pdf", pdf)
            out = mount_authorized_file(
                db_path=root / "library.sqlite",
                collection_id=cid,
                file_path=pdf,
                access_method="user_provided_file",
                access_origin="institutional",
                try_ocr=False,
            )
            self.assertTrue(out["ok"], out)
            self.assertEqual(out["access_method"], "user_provided_file")
            conn = open_db(root / "library.sqlite")
            row = conn.execute(
                "SELECT access_origin, access_method FROM pub_files WHERE publication_id=?",
                (out["publication_id"],),
            ).fetchone()
            conn.close()
            self.assertEqual(row["access_origin"], "institutional")
            self.assertEqual(row["access_method"], "user_provided_file")

    def test_fetch_mock_and_saml_guidance(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            out = institutional_fetch_and_mount(
                db_path=root / "library.sqlite",
                collection_id=cid,
                url="https://example.org/inst/paper.pdf",
                dest_dir=root / "papers",
                access_method="vpn_proxy",
                doi="10.1000/inst.mock.001",
                mock_fixture=FX / "tiny.pdf",
            )
            self.assertTrue(out["ok"], out)
            self.assertEqual(out["access_method"], "vpn_proxy")
            self.assertTrue(Path(out["saved_path"]).is_file())

            saml = institutional_fetch_and_mount(
                db_path=root / "library.sqlite",
                collection_id=cid,
                url="https://example.org/saml",
                dest_dir=root / "papers",
                access_method="saml_browser",
            )
            self.assertFalse(saml["ok"])
            self.assertEqual(saml.get("error"), "saml_browser_host_concern")

    def test_cli_requires_flag(self):
        script = ROOT / "skills" / "literature-library" / "scripts" / "institutional_access.py"
        proc = subprocess.run(
            [
                sys.executable,
                str(script),
                "--project_root",
                "/tmp",
                "--collection_id",
                "x",
                "--method",
                "user_provided_file",
                "--file_path",
                "/tmp/nope.pdf",
            ],
            capture_output=True,
            text=True,
        )
        # missing --institutional_access → argparse error
        self.assertNotEqual(proc.returncode, 0)


class McpClientTests(unittest.TestCase):
    def test_resolve_backend(self):
        self.assertEqual(resolve_backend("http"), "http")
        self.assertEqual(resolve_backend("mcp"), "mcp")
        old = os.environ.get("LIT_SEARCH_BACKEND")
        try:
            os.environ["LIT_SEARCH_BACKEND"] = "mcp"
            self.assertEqual(resolve_backend(None), "mcp")
        finally:
            if old is None:
                os.environ.pop("LIT_SEARCH_BACKEND", None)
            else:
                os.environ["LIT_SEARCH_BACKEND"] = old

    def test_search_and_related_mock(self):
        fx = FX / "mcp_mock_hits.json"
        hits, meta = search_publications_mcp(
            "exhaustion", sources=["pubmed"], limit=5, mock_fixture=fx
        )
        self.assertTrue(meta.get("ok"))
        self.assertEqual(meta.get("provider"), "mcp_mock")
        self.assertGreaterEqual(len(hits), 2)
        self.assertTrue(hits[0].get("doi", "").startswith("10.1000/mcp"))

        rhits, rmeta = related_publications_mcp(
            relation="similar", doi="10.1000/seed", mock_fixture=fx
        )
        self.assertTrue(rmeta.get("ok"))
        self.assertGreaterEqual(len(rhits), 1)

    def test_search_cli_mcp_mock(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cid = _col(root)
            # register run
            conn = open_db(root / "library.sqlite")
            rid = "run-mcp"
            rdir = root / "collections" / cid / "runs" / rid
            rdir.mkdir(parents=True)
            conn.execute(
                "INSERT INTO runs(id, collection_id, query, sources_json, started_at, status) "
                "VALUES (?,?,?,?,?,?)",
                (rid, cid, "mcp query", json.dumps(["pubmed"]), utc_now(), "configured"),
            )
            conn.commit()
            conn.close()
            (rdir / "run.json").write_text(
                json.dumps(
                    {
                        "query": "mcp query",
                        "sources": ["pubmed"],
                        "limit_per_source": 5,
                        "protocol_version": "0.6.0",
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "skills" / "literature-library" / "scripts" / "search.py"),
                    "--project_root",
                    str(root),
                    "--collection_id",
                    cid,
                    "--run_id",
                    rid,
                    "--backend",
                    "mcp",
                    "--mcp_mock_fixture",
                    str(FX / "mcp_mock_hits.json"),
                    "--skip_ingest",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            out = json.loads(proc.stdout.strip().splitlines()[-1])
            self.assertTrue(out["ok"])
            self.assertEqual(out.get("backend"), "mcp")
            self.assertGreaterEqual(out.get("n_hits", 0), 2)


class ResearchEvidenceBridgeTests(unittest.TestCase):
    def test_hits_jsonl_bridge(self):
        script = REEFV / "skills" / "literature-search" / "scripts" / "search.py"
        if not script.is_file():
            self.skipTest(
                "research-evidence is a separate plugin and is not beside this package"
            )
        with tempfile.TemporaryDirectory() as td:
            hits = Path(td) / "hits.jsonl"
            hits.write_text(
                json.dumps(
                    {
                        "doi": "10.1000/bridge.1",
                        "title": "Bridge hit",
                        "source": "mock",
                        "abstract": "x",
                        "url": "https://example.org/b1",
                    }
                )
                + "\n",
                encoding="utf-8",
            )
            log = Path(td) / "search_log.jsonl"
            proc = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--query",
                    "bridge",
                    "--hits_jsonl",
                    str(hits),
                    "--log_path",
                    str(log),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            out = json.loads(proc.stdout.strip().splitlines()[-1])
            self.assertTrue(out["ok"])
            self.assertTrue(out.get("delegated"))
            self.assertEqual(out.get("n_hits"), 1)
            self.assertTrue(log.is_file())

    def test_delegate_lit_mock_smoke(self):
        script = REEFV / "skills" / "literature-search" / "scripts" / "search.py"
        if not script.is_file():
            self.skipTest(
                "research-evidence is a separate plugin and is not beside this package"
            )
        proc = subprocess.run(
            [
                sys.executable,
                str(script),
                "--query",
                "exhaustion T cells",
                "--delegate_lit",
                "--lit_root",
                str(ROOT),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
        out = json.loads(proc.stdout.strip().splitlines()[-1])
        self.assertTrue(out["ok"], out)
        self.assertTrue(out.get("delegated"))
        self.assertGreaterEqual(out.get("n_hits", 0), 1)


if __name__ == "__main__":
    unittest.main()
