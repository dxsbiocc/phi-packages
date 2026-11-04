#!/usr/bin/env python3
"""Offline tests: HTTP version ties (mocked) + wiki source repoint (0.4.8)."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "tests" / "fixtures" / "providers"
sys.path.insert(0, str(ROOT / "skills" / "lit-db" / "scripts"))
sys.path.insert(0, str(ROOT / "skills" / "lit-search" / "scripts"))
sys.path.insert(0, str(ROOT / "skills" / "lit-wiki-lib" / "scripts"))

from db import (  # noqa: E402
    apply_version_detection,
    drain_id_remaps,
    init_db,
    link_to_collection,
    merge_version_proposals,
    open_db,
    proposals_from_doi_ties,
    resolve_or_create,
    utc_now,
)
import providers_http as ph  # noqa: E402
from wiki_lib import init_vault, safe_id, wiki_repoint_publication  # noqa: E402


class ExtractTiesTests(unittest.TestCase):
    def test_crossref_is_preprint_of(self):
        raw = json.loads((FX / "crossref_preprint_work.json").read_text())
        work = raw["message"]
        ties = ph.extract_crossref_version_ties(work)
        self.assertEqual(len(ties), 1)
        t = ties[0]
        self.assertEqual(t["from_doi"], "10.1101/2020.01.01.123456")
        self.assertEqual(t["to_doi"], "10.1038/s41586-020-0001-x")
        self.assertEqual(t["relation"], "preprint_of")
        self.assertTrue(t["high_confidence"])

    def test_openalex_locations_published_version(self):
        work = json.loads((FX / "openalex_preprint_work.json").read_text())
        ties = ph.extract_openalex_version_ties(work)
        self.assertGreaterEqual(len(ties), 1)
        t = ties[0]
        self.assertEqual(t["from_doi"], "10.1101/2020.01.01.123456")
        self.assertEqual(t["to_doi"], "10.1038/s41586-020-0001-x")
        self.assertTrue(t["high_confidence"])

    def test_discover_soft_fail(self):
        def boom(*a, **k):
            raise OSError("network down")

        with mock.patch.object(ph, "_get_json", side_effect=boom):
            ties, metas = ph.discover_version_ties(
                doi="10.1101/2020.01.01.123456", providers=["openalex", "crossref"]
            )
        self.assertEqual(ties, [])
        self.assertTrue(metas)
        self.assertTrue(all(m.get("provider") in ("error", "http") for m in metas))


class HttpProposalMappingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db = Path(self.tmp.name) / "library.sqlite"
        init_db(self.db)
        self.conn = open_db(self.db)

    def tearDown(self):
        self.conn.close()
        self.tmp.cleanup()

    def test_proposals_from_fixtures_and_apply(self):
        pre = resolve_or_create(
            self.conn,
            {
                "doi": "10.1101/2020.01.01.123456",
                "title": "Amazing Exhaustion Study In Mice And Humans",
                "authors": [{"family": "Smith"}],
                "year": 2020,
                "type": "preprint",
            },
        )
        jou = resolve_or_create(
            self.conn,
            {
                "doi": "10.1038/s41586-020-0001-x",
                "title": "Amazing Exhaustion Study In Mice And Humans",
                "authors": [{"family": "Smith"}],
                "year": 2020,
                "type": "journal-article",
            },
        )
        self.conn.commit()

        work = json.loads((FX / "openalex_preprint_work.json").read_text())
        ties = ph.extract_openalex_version_ties(work)
        props = proposals_from_doi_ties(self.conn, ties)
        self.assertGreaterEqual(len(props), 1)
        hit = props[0]
        self.assertEqual(hit["from_id"], pre["publication_id"])
        self.assertEqual(hit["to_id"], jou["publication_id"])
        self.assertEqual(hit["via"], "http")

        # Merge with offline (may also find soft match) then auto-merge
        out = apply_version_detection(
            self.conn, auto_merge_high_confidence=True, extra_proposals=props
        )
        self.conn.commit()
        self.assertTrue(out["ok"])
        self.assertGreaterEqual(out["n_merged"], 1)
        remaps = drain_id_remaps(self.conn)
        self.assertTrue(remaps)
        n = self.conn.execute("SELECT COUNT(*) c FROM publications").fetchone()["c"]
        self.assertEqual(n, 1)

    def test_merge_version_proposals_dedupes(self):
        a = [
            {
                "from_id": "doi:10.1101/a",
                "to_id": "doi:10.1/b",
                "relation": "preprint_of",
                "confidence": 0.7,
            }
        ]
        b = [
            {
                "from_id": "doi:10.1101/a",
                "to_id": "doi:10.1/b",
                "relation": "preprint_of",
                "confidence": 0.95,
                "high_confidence": True,
            }
        ]
        m = merge_version_proposals(a, b)
        self.assertEqual(len(m), 1)
        self.assertEqual(m[0]["confidence"], 0.95)


class WikiRepointTests(unittest.TestCase):
    def test_repoint_source_and_stub(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            conn = open_db(db)
            now = utc_now()
            col_id = "col-wiki"
            col_path = root / "collections" / col_id
            col_path.mkdir(parents=True)
            conn.execute(
                "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
                (col_id, "Wiki", "", now, str(col_path)),
            )
            old_id = "pmid:999888"
            new_id = "doi:10.1000/repoint.demo"
            # Create pubmed-only then upgrade path simulated
            r = resolve_or_create(
                conn,
                {"pmid": "999888", "title": "Repoint Demo Paper About Exhaustion", "year": 2021, "authors": ["A"]},
            )
            self.assertEqual(r["publication_id"], old_id)
            link_to_collection(conn, col_id, old_id)
            conn.commit()

            vault_info = init_vault(col_path, collection_id=col_id, collection_name="Wiki")
            self.assertTrue(vault_info["ok"])
            vault = col_path / "wiki-vault"
            old_safe = safe_id(old_id)
            new_safe = safe_id(new_id)
            sources = vault / "wiki" / "sources"
            sources.mkdir(parents=True, exist_ok=True)
            page = sources / f"{old_safe}.md"
            page.write_text(
                f"# Demo\n\n- publication_id: `{old_id}`\n- raw: [`publications/{old_safe}/body.txt`]"
                f"(../../raw/publications/{old_safe}/body.txt)\n",
                encoding="utf-8",
            )
            raw = vault / "raw" / "publications" / old_safe
            raw.mkdir(parents=True)
            (raw / "meta.json").write_text(
                json.dumps({"publication_id": old_id, "safe_id": old_safe}),
                encoding="utf-8",
            )
            (raw / "body.txt").write_text("body\n", encoding="utf-8")
            index = vault / "wiki" / "index.md"
            index.write_text(
                f"# Index\n\n- [Demo](sources/{old_safe}.md)\n",
                encoding="utf-8",
            )
            contrib = vault / "wiki" / "concepts" / "_contrib" / "exhaustion"
            contrib.mkdir(parents=True)
            (contrib / f"{old_safe}.json").write_text(
                json.dumps({"source_safe_id": old_safe, "publication_id": old_id, "note": "x"}),
                encoding="utf-8",
            )

            # Simulate DB reassign
            from db import _reassign_publication_id

            _reassign_publication_id(conn, old_id, new_id)
            conn.commit()
            remaps = drain_id_remaps(conn)
            self.assertEqual(remaps[0]["old_id"], old_id)
            self.assertEqual(remaps[0]["new_id"], new_id)

            out = wiki_repoint_publication(
                root, old_id, new_id, collection_id=col_id, db_path=db, leave_stub=True
            )
            self.assertTrue(out["ok"])
            self.assertGreaterEqual(out["n_vaults"], 1)
            new_page = sources / f"{new_safe}.md"
            self.assertTrue(new_page.is_file())
            body = new_page.read_text(encoding="utf-8")
            self.assertIn(f"publication_id: `{new_id}`", body)
            stub = page.read_text(encoding="utf-8")
            self.assertIn("Redirect", stub)
            self.assertIn(new_safe, stub)
            self.assertTrue((vault / "raw" / "publications" / new_safe / "meta.json").is_file())
            meta = json.loads(
                (vault / "raw" / "publications" / new_safe / "meta.json").read_text()
            )
            self.assertEqual(meta["publication_id"], new_id)
            self.assertFalse((contrib / f"{old_safe}.json").exists())
            self.assertTrue((contrib / f"{new_safe}.json").is_file())
            idx = index.read_text(encoding="utf-8")
            self.assertIn(f"sources/{new_safe}.md", idx)
            self.assertNotIn(f"sources/{old_safe}.md", idx)
            conn.close()


class CliMockHttpTests(unittest.TestCase):
    def test_versions_cli_mock_detect(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            conn = open_db(db)
            resolve_or_create(
                conn,
                {
                    "doi": "10.1101/2020.01.01.123456",
                    "title": "Amazing Exhaustion Study In Mice And Humans",
                    "year": 2020,
                    "authors": ["Smith"],
                    "type": "preprint",
                },
            )
            resolve_or_create(
                conn,
                {
                    "doi": "10.1038/s41586-020-0001-x",
                    "title": "Amazing Exhaustion Study In Mice And Humans",
                    "year": 2020,
                    "authors": ["Smith"],
                    "type": "journal-article",
                },
            )
            conn.commit()
            conn.close()

            out = subprocess.check_output(
                [
                    sys.executable,
                    str(ROOT / "skills/lit-versions/scripts/versions.py"),
                    "--project_root",
                    str(root),
                    "--action",
                    "detect",
                    "--mock",
                    "--fixture_dir",
                    str(FX),
                ],
                text=True,
            )
            data = json.loads(out)
            self.assertTrue(data["ok"])
            self.assertTrue(data["http_mock"])
            self.assertGreaterEqual(data["n_http"], 1)
            self.assertGreaterEqual(data["n"], 1)


if __name__ == "__main__":
    unittest.main()
