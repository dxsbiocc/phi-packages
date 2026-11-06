#!/usr/bin/env python3
"""Offline tests for lit-expand (mocked HTTP + confirm_write ingest)."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "literature-library" / "scripts"))

import providers_http as ph  # noqa: E402
import expand as expand_mod  # noqa: E402
from db import init_db, open_db  # noqa: E402

FX = ROOT / "tests" / "fixtures" / "providers"
EXPAND_PY = ROOT / "skills" / "literature-library" / "scripts" / "expand.py"
COLLECTION_PY = ROOT / "skills" / "literature-library" / "scripts" / "collection.py"


class FakeHTTPResponse:
    def __init__(self, data: bytes, status: int = 200):
        self._data = data
        self.status = status

    def read(self):
        return self._data

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def _patch_get(mapping: dict[str, bytes]):
    def fake_get(url: str, *, headers=None):  # noqa: ARG001
        for key, payload in mapping.items():
            if key in url:
                return payload
        raise AssertionError(f"unexpected URL in offline test: {url}")

    return mock.patch.object(ph, "_get", side_effect=fake_get)


class DedupTests(unittest.TestCase):
    def test_dedup_prefers_doi(self):
        hits = [
            {"doi": "10.1000/a", "pmid": "1", "title": "A"},
            {"doi": "10.1000/a", "pmid": "2", "title": "A dup"},
            {"pmid": "3", "title": "B"},
        ]
        out = ph.dedup_hits(hits, limit=10)
        self.assertEqual(len(out), 2)
        self.assertEqual(out[0]["doi"], "10.1000/a")


class PubmedNeighborsTests(unittest.TestCase):
    def test_elink_fixture(self):
        payload = (FX / "pubmed_elink.json").read_bytes()
        with _patch_get({"elink.fcgi": payload}):
            ids, meta = ph.pubmed_neighbors("28622514")
        self.assertEqual(meta["via"], "elink")
        self.assertEqual(ids, ["30092201", "31028727"])  # seed excluded


def _url_plain(url: str) -> str:
    from urllib.parse import unquote
    return unquote(url)


class ExpandRelatedTests(unittest.TestCase):
    def test_similar_via_pmid(self):
        elink = (FX / "pubmed_elink.json").read_bytes()
        efetch = (FX / "pubmed_efetch_neighbors.xml").read_bytes()
        with _patch_get({"elink.fcgi": elink, "efetch.fcgi": efetch}):
            hits, meta = ph.expand_related(relation="similar", limit=5, pmid="28622514")
        self.assertTrue(meta.get("ok"))
        self.assertEqual(meta["provider"], "pubmed")
        self.assertEqual(len(hits), 2)
        self.assertIn("Related T cell", hits[0]["title"])

    def test_references_via_openalex(self):
        detail = (FX / "openalex_work_detail.json").read_bytes()
        hydrate = (FX / "openalex_hydrate_refs.json").read_bytes()

        def fake_get(url: str, *, headers=None):  # noqa: ARG001
            u = _url_plain(url)
            if "filter=openalex:" in u or "openalex:W333" in u:
                return hydrate
            if "/works/W2741809807" in u or "/works/doi:" in u:
                return detail
            raise AssertionError(url)

        with mock.patch.object(ph, "_get", side_effect=fake_get):
            hits, meta = ph.expand_related(
                relation="references",
                limit=3,
                doi="10.1016/j.cell.2017.05.035",
                openalex="W2741809807",
            )
        self.assertTrue(meta.get("ok"))
        self.assertEqual(meta["provider"], "openalex")
        self.assertEqual(meta["via"], "referenced_works")
        self.assertEqual(len(hits), 3)
        self.assertTrue(any("Referenced work" in (h.get("title") or "") for h in hits))

    def test_citations_via_cites_filter(self):
        detail = (FX / "openalex_work_detail.json").read_bytes()
        cites = (FX / "openalex_cites.json").read_bytes()

        def fake_get(url: str, *, headers=None):  # noqa: ARG001
            u = _url_plain(url)
            if "cites:" in u:
                return cites
            if "/works/W2741809807" in u or "/works/doi:" in u:
                return detail
            raise AssertionError(url)

        with mock.patch.object(ph, "_get", side_effect=fake_get):
            hits, meta = ph.expand_related(
                relation="citations",
                limit=2,
                openalex="W2741809807",
                doi="10.1016/j.cell.2017.05.035",
            )
        self.assertTrue(meta.get("ok"))
        self.assertEqual(meta["via"], "cites_filter")
        self.assertEqual(len(hits), 2)
        self.assertIn("citing Zheng", hits[0]["title"])

    def test_similar_openalex_fallback_without_pmid(self):
        detail = (FX / "openalex_work_detail.json").read_bytes()
        # strip pmid from detail for this path — force OpenAlex related_works
        detail_obj = json.loads(detail)
        detail_obj["ids"].pop("pmid", None)
        detail = json.dumps(detail_obj).encode()
        hydrate = (FX / "openalex_hydrate_related.json").read_bytes()

        def fake_get(url: str, *, headers=None):  # noqa: ARG001
            u = _url_plain(url)
            if "filter=openalex:" in u or "openalex:W111" in u:
                return hydrate
            if "/works/W2741809807" in u or "/works/doi:" in u:
                return detail
            raise AssertionError(url)

        with mock.patch.object(ph, "_get", side_effect=fake_get):
            hits, meta = ph.expand_related(
                relation="similar",
                limit=2,
                doi="10.1016/j.cell.2017.05.035",
                openalex="W2741809807",
            )
        self.assertEqual(meta["via"], "related_works")
        self.assertEqual(len(hits), 2)


class ExpandCliTests(unittest.TestCase):
    def _make_project(self, tmp: Path) -> str:
        out = subprocess.check_output(
            [
                sys.executable,
                str(COLLECTION_PY),
                "--project_root",
                str(tmp),
                "--action",
                "create",
                "--name",
                "expand-test",
                "--no_wiki",
            ],
            text=True,
        )
        data = json.loads(out)
        self.assertTrue(data.get("ok"), data)
        return data["collection_id"]

    def test_mock_preview_no_write(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            cid = self._make_project(tmp)
            out = subprocess.check_output(
                [
                    sys.executable,
                    str(EXPAND_PY),
                    "--project_root",
                    str(tmp),
                    "--collection_id",
                    cid,
                    "--seed_id",
                    "28622514",
                    "--relation",
                    "similar",
                    "--limit",
                    "2",
                    "--mock",
                ],
                text=True,
            )
            data = json.loads(out)
            self.assertTrue(data["ok"])
            self.assertFalse(data.get("stub"))
            self.assertFalse(data["wrote"])
            self.assertEqual(data["n_hits"], 2)
            self.assertIsNone(data.get("hits_path"))
            # no expands dir
            cols = list((tmp / "collections").glob("*/expands"))
            self.assertEqual(cols, [])

    def test_mock_confirm_write_ingests(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            cid = self._make_project(tmp)
            out = subprocess.check_output(
                [
                    sys.executable,
                    str(EXPAND_PY),
                    "--project_root",
                    str(tmp),
                    "--collection_id",
                    cid,
                    "--seed_id",
                    "10.1016/j.cell.2017.05.035",
                    "--relation",
                    "references",
                    "--limit",
                    "2",
                    "--mock",
                    "--confirm_write",
                    "--no_enrich",
                ],
                text=True,
            )
            data = json.loads(out)
            self.assertTrue(data["ok"])
            self.assertTrue(data["wrote"])
            self.assertTrue(data.get("hits_path"))
            self.assertTrue(Path(data["hits_path"]).is_file())
            self.assertIsNotNone(data.get("ingest"))
            conn = open_db(tmp / "library.sqlite")
            try:
                n = conn.execute("SELECT COUNT(*) AS c FROM publications").fetchone()["c"]
            finally:
                conn.close()
            self.assertGreaterEqual(n, 2)

    def test_dry_run_blocks_confirm(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            cid = self._make_project(tmp)
            out = subprocess.check_output(
                [
                    sys.executable,
                    str(EXPAND_PY),
                    "--project_root",
                    str(tmp),
                    "--collection_id",
                    cid,
                    "--seed_id",
                    "28622514",
                    "--mock",
                    "--confirm_write",
                    "--dry-run",
                ],
                text=True,
            )
            data = json.loads(out)
            self.assertFalse(data["wrote"])

    def test_parse_seed_variants(self):
        self.assertEqual(expand_mod._parse_seed_string("PMID:28622514")["pmid"], "28622514")
        self.assertEqual(
            expand_mod._parse_seed_string("10.1016/j.cell.2017.05.035")["doi"],
            "10.1016/j.cell.2017.05.035",
        )
        self.assertEqual(expand_mod._parse_seed_string("W2741809807")["openalex"], "W2741809807")
        self.assertEqual(expand_mod._parse_seed_string("PMC5533227")["pmcid"], "PMC5533227")


if __name__ == "__main__":
    unittest.main()
