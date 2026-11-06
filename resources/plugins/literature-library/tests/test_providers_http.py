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
sys.path.insert(0, str(ROOT / "skills" / "lit-search" / "scripts"))

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
        self.assertEqual(ph.normalize_source("S2"), "semantic_scholar")
        self.assertEqual(ph.normalize_source("semantic_scholar"), "semantic_scholar")
        self.assertEqual(ph.normalize_source("epmc"), "europepmc")
        self.assertEqual(ph.normalize_source("pubmed_central"), "pmc")
        self.assertEqual(ph.normalize_source("data-cite"), "datacite")


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


class SemanticScholarTests(unittest.TestCase):
    def test_s2_alias_and_parse(self):
        payload = (FX / "s2_search.json").read_bytes()
        with _patch_get({"api.semanticscholar.org/graph/v1/paper/search": payload}):
            hits, meta = ph.search_source("S2", "single-cell", limit=5)
        self.assertEqual(meta["provider"], "http")
        self.assertEqual(meta["source"], "semantic_scholar")
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["title"], "Semantic Scholar fixture paper")
        self.assertEqual(hits[0]["doi"], "10.1000/s2.fixture")

    def test_s2_soft_fail_401(self):
        def boom(url, *, headers=None):  # noqa: ARG001
            raise urllib.error.HTTPError(url, 401, "Unauthorized", hdrs=None, fp=None)

        with mock.patch.object(ph, "_get", side_effect=boom):
            hits, meta = ph.search_source("semantic_scholar", "anything", limit=3)
        self.assertEqual(hits, [])
        self.assertEqual(meta["provider"], "error")
        self.assertIn("S2_API_KEY", meta["message"])


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



class PmcTests(unittest.TestCase):
    def test_search_pmc_fixture(self):
        esearch = (FX / "pmc_esearch.json").read_bytes()
        esummary = (FX / "pmc_esummary.json").read_bytes()
        with _patch_get(
            {
                "eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi": esearch,
                "eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi": esummary,
            }
        ):
            hits, meta = ph.search_source("pmc", "infiltrating T cells liver", limit=5)
        self.assertEqual(meta["provider"], "http")
        self.assertEqual(meta["source"], "pmc")
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["pmcid"], "PMC5508784")
        self.assertEqual(hits[0]["pmid"], "28622514")
        self.assertEqual(hits[0]["doi"], "10.1016/j.cell.2017.05.035")
        self.assertEqual(hits[0]["source"], "pmc")
        self.assertIn("Infiltrating T Cells", hits[0]["title"])

    def test_hit_from_pmc_unit(self):
        data = json.loads((FX / "pmc_esummary.json").read_text())
        row = data["result"]["5508784"]
        h = ph._hit_from_pmc_esummary("5508784", row)
        self.assertEqual(h["pmcid"], "PMC5508784")
        self.assertEqual(h["year"], "2017")


class DoajTests(unittest.TestCase):
    def test_search_doaj_fixture(self):
        payload = (FX / "doaj_articles.json").read_bytes()
        with _patch_get({"doaj.org/api/search/articles": payload}):
            hits, meta = ph.search_source("doaj", "infiltrating", limit=5)
        self.assertEqual(meta["provider"], "http")
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["source"], "doaj")
        self.assertTrue(hits[0]["doi"])
        self.assertTrue(hits[0]["title"])
        self.assertEqual(hits[0]["type"], "journal-article")

    def test_doaj_soft_fail_shape(self):
        def boom(url, *, headers=None):  # noqa: ARG001
            raise urllib.error.HTTPError(url, 500, "Server Error", hdrs=None, fp=None)

        with mock.patch.object(ph, "_get", side_effect=boom):
            hits, meta = ph.search_source("doaj", "x", limit=1)
        self.assertEqual(hits, [])
        self.assertEqual(meta["provider"], "error")
        self.assertEqual(meta["source"], "doaj")


class DataciteTests(unittest.TestCase):
    def test_search_datacite_fixture(self):
        payload = (FX / "datacite_dois.json").read_bytes()
        with _patch_get({"api.datacite.org/dois": payload}):
            hits, meta = ph.search_source("datacite", "hepatocellular", limit=5)
        self.assertEqual(meta["provider"], "http")
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["source"], "datacite")
        self.assertTrue(hits[0]["doi"])
        self.assertTrue(hits[0]["title"])

    def test_datacite_unit(self):
        data = json.loads((FX / "datacite_dois.json").read_text())
        h = ph._hit_from_datacite(data["data"][0])
        self.assertIsNotNone(h)
        self.assertEqual(h["source"], "datacite")


class ZenodoTests(unittest.TestCase):
    def test_search_zenodo_fixture(self):
        payload = (FX / "zenodo_records.json").read_bytes()
        with _patch_get({"zenodo.org/api/records": payload}):
            hits, meta = ph.search_source("zenodo", "hepatocellular", limit=5)
        self.assertEqual(meta["provider"], "http")
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["source"], "zenodo")
        self.assertIn("10.5281/zenodo", hits[0]["doi"] or "")
        self.assertIn("Hepatocellular", hits[0]["title"])

    def test_zenodo_soft_fail(self):
        def boom(url, *, headers=None):  # noqa: ARG001
            raise urllib.error.URLError("timed out")

        with mock.patch.object(ph, "_get", side_effect=boom):
            hits, meta = ph.search_source("zenodo", "x", limit=1)
        self.assertEqual(hits, [])
        self.assertEqual(meta["provider"], "error")



if __name__ == "__main__":
    unittest.main()
