"""Regression tests for problems found by running the plugin on real PDFs and live APIs."""
from __future__ import annotations

import os
import unittest
from unittest import mock

import ingest
import oa_verify
import providers_http
from enrich import needs_enrich
from pdf_extract import guess_abstract, pdf_info_title


class PdfTextTests(unittest.TestCase):
    def test_wrapped_doi_is_rejoined(self):
        ids = ingest.extract_ids_from_text("https://doi.org/10.1038/\ns41575-020-0342-4\n\nThe unique")
        self.assertEqual(ids["doi"], "10.1038/s41575-020-0342-4")

    def test_abstract_is_the_paragraph_not_the_body(self):
        body = "Title\n\nAbstract | " + ("Liver regeneration restores mass. " * 12) + "\n\nIntroduction\nBody text."
        got = guess_abstract(body)
        self.assertTrue(got.startswith("Liver regeneration"))
        self.assertNotIn("Body text", got)
        self.assertIsNone(guess_abstract("no heading here " * 50))

    def test_pdf_metadata_title_rejects_site_names(self):
        with mock.patch("pdf_extract.subprocess.run") as run, mock.patch("pdf_extract.shutil.which", return_value="/x/pdfinfo"):
            run.return_value.stdout = b"Title:  Science Journals \xe2\x80\x94 AAAS\nPages: 3\n"
            self.assertIsNone(pdf_info_title(__import__("pathlib").Path("a.pdf")))
            run.return_value.stdout = b"Title:  An Agentic AI Framework for Single-Cell Data\n"
            self.assertEqual(pdf_info_title(__import__("pathlib").Path("a.pdf")), "An Agentic AI Framework for Single-Cell Data")


class EnrichGateTests(unittest.TestCase):
    def test_paper_with_doi_but_no_abstract_is_enriched(self):
        self.assertTrue(needs_enrich({"title": "A real title here", "doi": "10.1/x", "abstract": ""}))

    def test_already_looked_up_paper_is_not_retried(self):
        row = {"title": "A real title here", "doi": "10.1/x", "abstract": "", "enriched_at": "2026-01-01T00:00:00Z"}
        self.assertFalse(needs_enrich(row))

    def test_paper_with_abstract_is_left_alone(self):
        row = {"title": "T title", "doi": "10.1/x", "abstract": "abs", "authors_json": '["A, B"]'}
        self.assertFalse(needs_enrich(row))

    def test_pdf_record_without_authors_is_enriched_even_with_an_abstract(self):
        self.assertTrue(needs_enrich({"title": "GraphGPS", "arxiv": "2205.12454", "abstract": "guessed", "authors_json": None}))

    def test_no_id_but_a_real_title_is_matched_by_title(self):
        self.assertTrue(needs_enrich({"title": "Causal machine learning for single-cell genomics", "abstract": ""}))

    def test_no_id_and_a_short_title_is_not_looked_up(self):
        self.assertFalse(needs_enrich({"title": "GraphGPS", "abstract": ""}))

    def test_arxiv_only_paper_is_enriched(self):
        self.assertTrue(needs_enrich({"title": "GraphGPS", "arxiv": "2205.12454", "abstract": ""}))


class OaGateTests(unittest.TestCase):
    def test_live_classification_downgrades_html_on_pdf_url(self):
        class Resp:
            status = 200
            headers = {"Content-Type": "text/html; charset=utf-8"}
            def getcode(self): return 200
            def geturl(self): return "https://x.org/a.pdf"
            def __enter__(self): return self
            def __exit__(self, *a): return False
        out = oa_verify.verify_oa_url("https://x.org/a.pdf", opener=lambda req, timeout: Resp())
        self.assertEqual(out["oa_status"], "unknown")
        self.assertIn("browser", out["message"])

    def test_html_landing_page_is_still_not_pdf(self):
        class Resp:
            status = 200
            headers = {"Content-Type": "text/html"}
            def getcode(self): return 200
            def geturl(self): return "https://x.org/article/123"
            def __enter__(self): return self
            def __exit__(self, *a): return False
        out = oa_verify.verify_oa_url("https://x.org/article/123", opener=lambda req, timeout: Resp())
        self.assertEqual(out["oa_status"], "not_pdf")


class TitleMatchTests(unittest.TestCase):
    def test_near_identical_title_is_accepted_and_different_one_rejected(self):
        item = {"DOI": "10.1038/s41588-025-02124-2", "title": ["Causal machine learning for single-cell genomics"], "issued": {"date-parts": [[2025]]}}
        with mock.patch.object(providers_http, "_get_json", return_value={"message": {"items": [item]}}):
            hit, meta = providers_http.find_by_title("Causal Machine Learning for Single-Cell Genomics.")
            self.assertEqual(hit["doi"], "10.1038/s41588-025-02124-2")
            hit, meta = providers_http.find_by_title("A completely different paper about liver regeneration")
            self.assertIsNone(hit)

    def test_short_title_is_not_searched(self):
        with mock.patch.object(providers_http, "_get_json") as get:
            hit, _ = providers_http.find_by_title("GraphGPS")
        self.assertIsNone(hit)
        get.assert_not_called()


class HttpTests(unittest.TestCase):
    def test_ncbi_urls_are_identified_and_keyed(self):
        with mock.patch.dict(os.environ, {"NCBI_API_KEY": "k123", "CONTACT_EMAIL": "me@example.org"}, clear=False):
            url = providers_http._ncbi_url("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed")
        self.assertIn("tool=phi-literature-library", url)
        self.assertIn("api_key=k123", url)
        self.assertIn("email=me%40example.org", url)

    def test_pubmed_search_sorts_by_relevance(self):
        seen = []
        with mock.patch.object(providers_http, "_get", side_effect=lambda u, **k: (seen.append(u), b'{"esearchresult":{"idlist":[]}}')[1]):
            providers_http.search_pubmed("liver", limit=3)
        self.assertNotIn("sort=", seen[0])

    def test_opencitations_uses_the_current_api(self):
        seen = []
        with mock.patch.object(providers_http, "_get_json", side_effect=lambda u, **k: (seen.append(u), [])[1]):
            providers_http.opencitations_dois("10.1126/science.aaf2403", "citations")
        self.assertEqual(seen[0], "https://api.opencitations.net/index/v1/citations/10.1126/science.aaf2403")

    def test_doi_lookup_falls_back_to_crossref(self):
        work = {"DOI": "10.1016/j.x.1", "title": ["Cloning of a gene"], "issued": {"date-parts": [[2020]]},
                "container-title": ["Aquaculture Reports"], "type": "journal-article"}
        with mock.patch.object(providers_http, "search_source", return_value=([], {"provider": "http"})), \
                mock.patch.object(providers_http, "get_crossref_work", return_value=(work, {"provider": "http", "source": "crossref"})):
            hit, meta = providers_http.fetch_by_doi("10.1016/j.x.1")
        self.assertEqual(meta["via"], "doi_crossref")
        self.assertEqual(hit["title"], "Cloning of a gene")

    def test_filename_like_title_is_a_placeholder(self):
        self.assertTrue(needs_enrich({"title": "Chou_2007_Spot14_receptor", "doi": "10.1/x", "abstract": "has one"}))

    def test_verify_oa_looks_up_openalex_when_record_has_no_link(self):
        with mock.patch.object(providers_http, "openalex_oa_urls", return_value=["https://x.org/a.pdf"]) as look, \
                mock.patch.object(oa_verify, "verify_oa_url", return_value={"ok": True, "oa_status": "pdf", "url": "https://x.org/a.pdf"}):
            out = oa_verify.verify_record_oa({"doi": "10.1/x"})
        look.assert_called_once()
        self.assertEqual(out["oa_status"], "pdf")


if __name__ == "__main__":
    unittest.main()
