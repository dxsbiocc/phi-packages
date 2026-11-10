#!/usr/bin/env python3
"""Unit tests for reference formats and the wiki payload validator."""
from __future__ import annotations

import json
import unittest

from formats import cite_keys, latex_to_unicode, parse_authors, to_bibtex, to_csl_json, to_ris
from wiki import slugify, validate_compile_payload

PUBS = [
    {"id": "doi:10.1/a", "title": "T cell exhaustion & the tumour_microenvironment", "year": 2020,
     "authors_json": json.dumps(["Zheng, Chunhong", "Jane Q Public"]), "venue": "Cell", "type": "journal-article",
     "doi": "10.1/a", "pmid": "1", "abstract": "100% of cells"},
    {"id": "doi:10.1/b", "title": "T cell exhaustion revisited", "year": 2020,
     "authors_json": json.dumps(["Zheng, Chunhong"]), "venue": "", "type": "", "doi": "10.1/b"},
    {"id": "arxiv:2001.1", "title": "Ünïcode études", "year": None, "authors_json": None, "arxiv": "2001.1"},
]


class FormatTests(unittest.TestCase):
    def test_authors_in_both_orders(self):
        a = parse_authors(json.dumps(["Zheng, Chunhong", "Jane Q Public", "Plato"]))
        self.assertEqual(a[0], {"family": "Zheng", "given": "Chunhong"})
        self.assertEqual(a[1], {"family": "Public", "given": "Jane Q"})
        self.assertEqual(a[2], {"family": "Plato", "given": ""})

    def test_europepmc_family_then_initials(self):
        a = parse_authors(["Luo H", "Grünewald TGP", "Gerstung M.", "Van Der Berg J"])
        self.assertEqual(a[0], {"family": "Luo", "given": "H"})
        self.assertEqual(a[1], {"family": "Grünewald", "given": "TGP"})
        self.assertEqual(a[2], {"family": "Gerstung", "given": "M"})
        self.assertEqual(a[3], {"family": "Van Der Berg", "given": "J"})

    def test_latex_accents_are_decoded(self):
        self.assertEqual(latex_to_unicode("St{\\aa}hl, Patrik L"), "Ståhl, Patrik L")
        self.assertEqual(latex_to_unicode('Gr\\"{u}newald'), "Grünewald")
        self.assertEqual(latex_to_unicode("Garc\\'{\\i}a"), "García")
        self.assertEqual(latex_to_unicode("{DNA} and {RNA} \\v{S}ebesta"), "DNA and RNA Šebesta")

    def test_cite_keys_are_unique_and_ascii(self):
        keys = cite_keys(PUBS)
        self.assertEqual(len(set(keys.values())), 3)
        self.assertEqual(keys["doi:10.1/a"], "zheng2020cell")
        self.assertEqual(keys["doi:10.1/b"], "zheng2020cella")  # same author, year and word: letter suffix
        self.assertTrue(all(k.isascii() for k in keys.values()))

    def test_bibtex_escapes_and_types(self):
        bib = to_bibtex(PUBS)
        self.assertIn(r"tumour\_microenvironment", bib)
        self.assertIn(r"100\% of cells", bib)
        self.assertIn("@article{zheng2020cell,", bib)
        self.assertIn("@misc{", bib)
        self.assertEqual(bib.count("@"), 3)

    def test_ris_and_csl(self):
        ris = to_ris(PUBS)
        self.assertEqual(ris.count("ER  - "), 3)
        self.assertIn("AU  - Zheng, Chunhong", ris)
        csl = json.loads(to_csl_json(PUBS))
        self.assertEqual(csl[0]["issued"], {"date-parts": [[2020]]})
        self.assertNotIn("issued", csl[2])


class WikiPayloadTests(unittest.TestCase):
    good = {"publication_id": "doi:1", "summary": "x" * 50, "key_claims": ["c"]}

    def test_accepts_minimal_and_folds_entities_into_concepts(self):
        out = validate_compile_payload({**self.good, "concepts": ["A"], "entities": [{"symbol": "LAYN"}, "a"]})
        self.assertEqual([c["name"] for c in out["concepts"]], ["A", "LAYN"])

    def test_rejects_short_summary_and_missing_claims(self):
        with self.assertRaises(ValueError):
            validate_compile_payload({**self.good, "summary": "short"})
        with self.assertRaises(ValueError):
            validate_compile_payload({**self.good, "key_claims": []})
        with self.assertRaises(ValueError):
            validate_compile_payload("not an object")

    def test_slugify_keeps_gene_names_whole(self):
        self.assertEqual(slugify("PD-L1"), "pd-l1")
        self.assertEqual(slugify("CD8+ T cells"), "cd8-plus-t-cells")


if __name__ == "__main__":
    unittest.main()
