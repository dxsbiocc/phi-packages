#!/usr/bin/env python3
"""Offline tests for lit-search HTTP provider adapters (mocked urllib)."""
from __future__ import annotations

import io
import json
import sys
import unittest
import urllib.error
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]

import providers_http as ph  # noqa: E402

FX = ROOT / "tests" / "fixtures" / "providers"


def _bytes_response(data: bytes):
    return io.BytesIO(data)


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
    """Patch ph._get to return fixture bytes keyed by URL substring."""

    def fake_get(url: str, *, headers=None):  # noqa: ARG001
        for key, payload in mapping.items():
            if key in url:
                return payload
        raise AssertionError(f"unexpected URL in offline test: {url}")

    return mock.patch.object(ph, "_get", side_effect=fake_get)


class ParseHelpersTests(unittest.TestCase):
    def test_inverted_abstract(self):
        idx = {"Hello": [0], "world": [1]}
        self.assertEqual(ph._inverted_abstract(idx), "Hello world")

    def test_doi_norm(self):
        self.assertEqual(ph._doi_norm("https://doi.org/10.1016/j.cell.2017.05.035"), "10.1016/j.cell.2017.05.035")
        self.assertIsNone(ph._doi_norm("not-a-doi"))

    def test_normalize_source_aliases(self):
        self.assertEqual(ph.normalize_source("epmc"), "europepmc")
        self.assertEqual(ph.normalize_source("Europe_PMC"), "europepmc")
        self.assertEqual(ph.normalize_source("pubmed"), "pubmed")
        self.assertNotIn("semantic_scholar", ph.SUPPORTED_SOURCES)


class OpenAlexTests(unittest.TestCase):
    def test_search_openalex_fixture(self):
        payload = (FX / "openalex_works.json").read_bytes()
        with _patch_get({"api.openalex.org/works": payload}):
            hits, meta = ph.search_source("openalex", "T cell liver cancer", limit=5)
        self.assertEqual(meta["provider"], "http")
        self.assertEqual(len(hits), 1)
        self.assertIn("Landscape of Infiltrating T Cells", hits[0]["title"])
        self.assertEqual(hits[0]["doi"], "10.1016/j.cell.2017.05.035")
        self.assertEqual(hits[0]["pmid"], "28622514")
        self.assertEqual(hits[0]["source"], "openalex")
        self.assertIn("T", hits[0]["abstract"] or "")


class CrossrefTests(unittest.TestCase):
    def test_search_crossref_fixture(self):
        payload = (FX / "crossref_works.json").read_bytes()
        with _patch_get({"api.crossref.org/works": payload}):
            hits, meta = ph.search_source("crossref", "HIV gp120", limit=5)
        self.assertEqual(meta["provider"], "http")
        self.assertEqual(len(hits), 1)
        self.assertIn("HIV-1", hits[0]["title"])
        self.assertEqual(hits[0]["doi"], "10.1038/nature12373")
        self.assertEqual(hits[0]["year"], "2013")
        self.assertNotIn("<p>", hits[0].get("abstract") or "")


class ArxivTests(unittest.TestCase):
    def test_search_arxiv_fixture(self):
        payload = (FX / "arxiv_feed.xml").read_bytes()
        with _patch_get({"export.arxiv.org/api/query": payload}):
            hits, meta = ph.search_source("arxiv", "attention transformer", limit=5)
        self.assertEqual(meta["provider"], "http")
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["title"], "Attention Is All You Need")
        self.assertEqual(hits[0]["arxiv"], "1706.03762")
        self.assertEqual(hits[0]["type"], "preprint")
        self.assertEqual(hits[0]["venue"], "arXiv")


class BiorxivTests(unittest.TestCase):
    def test_biorxiv_doi_detail(self):
        payload = (FX / "biorxiv_details.json").read_bytes()
        with _patch_get({"api.biorxiv.org/details/biorxiv/10.1101": payload}):
            hits, meta = ph.search_source(
                "biorxiv", "10.1101/2020.01.01.123456", limit=5
            )
        self.assertEqual(meta.get("via"), "doi_detail")
        self.assertEqual(len(hits), 1)
        self.assertIn("T cell exhaustion", hits[0]["title"])
        self.assertEqual(hits[0]["source"], "biorxiv")
        self.assertEqual(hits[0]["type"], "preprint")

    def test_biorxiv_keyword_via_epmc(self):
        payload = (FX / "epmc_biorxiv.json").read_bytes()
        with _patch_get({"europepmc/webservices/rest/search": payload}):
            hits, meta = ph.search_source("biorxiv", "T cell exhaustion", limit=5)
        self.assertEqual(meta.get("via"), "europepmc_journal")
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["source"], "biorxiv")
        self.assertIn("EPMC bioRxiv", hits[0]["title"])




class DispatcherUnknownTests(unittest.TestCase):
    def test_unknown_source(self):
        hits, meta = ph.search_source("google_scholar", "x", limit=1)
        self.assertEqual(hits, [])
        self.assertEqual(meta["provider"], "error")


class HitFromParsersUnit(unittest.TestCase):
    def test_openalex_unit(self):
        data = json.loads((FX / "openalex_works.json").read_text())
        h = ph._hit_from_openalex(data["results"][0])
        self.assertIsNotNone(h)
        self.assertEqual(h["openalex"], "W2741809807")

    def test_crossref_unit(self):
        data = json.loads((FX / "crossref_works.json").read_text())
        h = ph._hit_from_crossref(data["message"]["items"][0])
        self.assertEqual(h["venue"], "Nature")

    def test_biorxiv_unit(self):
        data = json.loads((FX / "biorxiv_details.json").read_text())
        h = ph._hit_from_biorxiv(data["collection"][0], "biorxiv")
        self.assertEqual(h["doi"], "10.1101/2020.01.01.123456")












if __name__ == "__main__":
    unittest.main()
