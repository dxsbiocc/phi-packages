#!/usr/bin/env python3
"""Offline UX polish tests for 0.5.3-draft."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

from db import (  # noqa: E402
    init_db,
    open_db,
    record_pub_file,
    utc_now,
)
from rebind_files import rebind_pub_files  # noqa: E402
from wiki import init_vault  # noqa: E402


FIX = ROOT / "tests" / "fixtures"




class QueryAliasTests(unittest.TestCase):
    def test_query_alias_smoke(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            cid = "col-q"
            col = root / "collections" / cid
            col.mkdir(parents=True)
            conn = open_db(db)
            conn.execute(
                "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
                (cid, "Q", "", utc_now(), str(col)),
            )
            conn.execute(
                "INSERT INTO publications(id,title,abstract,created_at,updated_at) VALUES (?,?,?,?,?)",
                ("doi:10.1/x", "Exhaustion study", "T cell exhaustion LAYN", utc_now(), utc_now()),
            )
            conn.execute(
                "INSERT INTO collection_pubs(collection_id,publication_id,added_at) VALUES (?,?,?)",
                (cid, "doi:10.1/x", utc_now()),
            )
            conn.commit()
            conn.close()
            init_vault(col, collection_id=cid, collection_name="Q")

            q_script = ROOT / "skills" / "manage-library" / "scripts" / "query.py"
            for flag in ("--query", "--q"):
                proc = subprocess.run(
                    [sys.executable, str(q_script), "--project_root", str(root),
                     "--collection_id", cid, flag, "exhaustion"],
                    capture_output=True, text=True,
                )
                self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
                out = json.loads(proc.stdout)
                self.assertTrue(out.get("ok"))
                self.assertIn("papers", out)


class RebindFilesTests(unittest.TestCase):
    def test_rebind_dry_and_confirm(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            old_dir = root / "old_papers"
            new_dir = root / "new_papers"
            old_dir.mkdir()
            new_dir.mkdir()
            old_path = old_dir / "paper.pdf"
            payload = b"%PDF-1.4 offline-rebind-fixture-bytes"
            old_path.write_bytes(payload)
            # hash while file exists
            import hashlib

            digest = hashlib.sha256(payload).hexdigest()
            conn = open_db(db)
            conn.execute(
                "INSERT INTO publications(id,title,created_at,updated_at) VALUES (?,?,?,?)",
                ("doi:10.rebind/1", "Rebind paper", utc_now(), utc_now()),
            )
            record_pub_file(
                conn,
                publication_id="doi:10.rebind/1",
                file_path=str(old_path),
                content_hash=digest,
                ingest_origin="local_file",
            )
            conn.commit()
            # move file
            new_path = new_dir / "paper_moved.pdf"
            old_path.rename(new_path)
            self.assertFalse(old_path.exists())

            dry = rebind_pub_files(conn, scan_dir=new_dir, confirm=False)
            self.assertTrue(dry["ok"])
            self.assertTrue(dry["dry_run"])
            self.assertEqual(dry["n_matched"], 1)
            self.assertEqual(dry["n_updated"], 0)
            row = conn.execute("SELECT file_path FROM pub_files").fetchone()
            self.assertEqual(row["file_path"], str(old_path))

            done = rebind_pub_files(conn, scan_dir=new_dir, confirm=True)
            self.assertEqual(done["n_updated"], 1)
            row = conn.execute("SELECT file_path, content_hash FROM pub_files").fetchone()
            self.assertEqual(row["file_path"], str(new_path.resolve()))
            self.assertEqual(row["content_hash"], digest)
            conn.close()








if __name__ == "__main__":
    unittest.main()
