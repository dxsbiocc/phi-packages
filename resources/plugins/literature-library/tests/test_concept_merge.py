#!/usr/bin/env python3
"""Offline tests: per-source concept contribs + mock deep merge (v0.4.5)."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "lit-wiki-lib" / "scripts"))
sys.path.insert(0, str(ROOT / "skills" / "lit-db" / "scripts"))

from wiki_lib import (  # noqa: E402
    ingest_publication_to_wiki,
    init_vault,
    list_concept_contributions,
    slugify,
)
from wiki_llm import (  # noqa: E402
    apply_concept_merge,
    build_concept_merge_request,
    mock_concept_merge_payload,
    validate_concept_merge_payload,
)


FIX = ROOT / "tests" / "fixtures"
ZHENG_ID = "doi:10.1016/j.cell.2017.05.035"
REL_ID = "doi:10.1038/s41586-example-exhaustion"


def _zheng_pub() -> dict:
    return {
        "id": ZHENG_ID,
        "title": "Landscape of Infiltrating T Cells in Liver Cancer Revealed by Single-Cell Sequencing",
        "abstract": (
            "Here we perform deep single-cell RNA sequencing on 5,063 single T cells from six HCC patients. "
            "We identify 11 T cell subsets. Exhausted CD8+ T cells and Tregs are enriched. "
            "LAYN (layilin) is upregulated on immunosuppressive T cells."
        ),
        "year": 2017,
        "doi": "10.1016/j.cell.2017.05.035",
        "pmid": "28622514",
        "needs_review": 0,
        "tags_json": "[]",
    }


def _related_pub() -> dict:
    return {
        "id": REL_ID,
        "title": "Shared programs of T cell exhaustion across human tumors (fixture)",
        "abstract": (
            "Exhausted CD8+ T cells accumulate in tumors with progressive dysfunction. "
            "Inhibitory receptors mark late exhausted states relevant to checkpoint blockade."
        ),
        "year": 2019,
        "doi": "10.1038/s41586-example-exhaustion",
        "needs_review": 0,
        "tags_json": "[]",
    }


class ConceptMergeTests(unittest.TestCase):
    def test_second_paper_writes_contrib_and_pending(self):
        with tempfile.TemporaryDirectory() as td:
            col = Path(td) / "collections" / "c1"
            init_vault(col, collection_id="c1", collection_name="merge-test")
            vault = col / "wiki-vault"
            r1 = ingest_publication_to_wiki(
                vault, _zheng_pub(), compile_mode="mock", llm_json=FIX / "mock_compile_zheng.json"
            )
            self.assertTrue(r1["ok"])
            self.assertEqual(r1.get("concept_merge_pending") or [], [])
            slug = slugify("T cell exhaustion")
            contribs1 = list_concept_contributions(vault, slug)
            self.assertEqual(len(contribs1), 1)
            page = vault / "wiki" / "concepts" / f"{slug}.md"
            before = page.read_text(encoding="utf-8")
            self.assertIn("## Definition", before)
            self.assertIn("## Sources", before)

            r2 = ingest_publication_to_wiki(
                vault,
                _related_pub(),
                compile_mode="mock",
                llm_json=FIX / "mock_compile_related_exhaustion.json",
            )
            self.assertTrue(r2["ok"])
            self.assertIn(slug, r2.get("concept_merge_pending") or [])
            contribs2 = list_concept_contributions(vault, slug)
            self.assertEqual(len(contribs2), 2)
            stub = vault / "raw" / f"concept_merge_request_{slug}.json"
            self.assertTrue(stub.is_file())
            stub_obj = json.loads(stub.read_text(encoding="utf-8"))
            self.assertEqual(stub_obj["action"], "concept_merge_pending")
            after_upsert = page.read_text(encoding="utf-8")
            # Still upsert-link Sources (both papers) before merge apply
            self.assertGreaterEqual(after_upsert.count("- ["), 2)

    def test_mock_merge_updates_definition_both_sources(self):
        with tempfile.TemporaryDirectory() as td:
            col = Path(td) / "collections" / "c1"
            init_vault(col, collection_id="c1")
            vault = col / "wiki-vault"
            ingest_publication_to_wiki(
                vault, _zheng_pub(), compile_mode="mock", llm_json=FIX / "mock_compile_zheng.json"
            )
            ingest_publication_to_wiki(
                vault,
                _related_pub(),
                compile_mode="mock",
                llm_json=FIX / "mock_compile_related_exhaustion.json",
            )
            slug = slugify("T cell exhaustion")
            req = build_concept_merge_request(vault, slug, response_path="/tmp/x.json")
            self.assertEqual(req["action"], "concept_merge_request")
            self.assertGreaterEqual(req["n_contributions"], 2)

            payload = mock_concept_merge_payload(vault, slug)
            self.assertIn("definition", payload)
            # definition should mention both source titles or ids somehow
            defn = payload["definition"]
            self.assertTrue(
                "Zheng" in defn
                or "Landscape" in defn
                or "Infiltrating" in defn
                or "fixture" in defn.lower()
                or "Shared programs" in defn
                or "2 compiled" in defn
            )
            result = apply_concept_merge(vault, slug, payload, backend="mock")
            self.assertTrue(result["ok"])
            md = (vault / "wiki" / "concepts" / f"{slug}.md").read_text(encoding="utf-8")
            self.assertIn("## Across papers", md)
            self.assertIn("## Definition", md)
            self.assertGreaterEqual(md.count("- ["), 2)  # Sources preserved
            self.assertIn("merge:", md)
            # synthesis bullets cite sources
            self.assertIn("`", md)
            log = (vault / "wiki" / "log.md").read_text(encoding="utf-8")
            self.assertIn("concept_merge", log)

    def test_needs_review_requires_banner(self):
        with tempfile.TemporaryDirectory() as td:
            col = Path(td) / "collections" / "c1"
            init_vault(col, collection_id="c1")
            vault = col / "wiki-vault"
            p1 = _zheng_pub()
            p2 = _related_pub()
            p2["needs_review"] = 1
            ingest_publication_to_wiki(
                vault, p1, compile_mode="mock", llm_json=FIX / "mock_compile_zheng.json"
            )
            ingest_publication_to_wiki(
                vault, p2, compile_mode="mock", llm_json=FIX / "mock_compile_related_exhaustion.json"
            )
            slug = slugify("T cell exhaustion")
            # mock path always sets banner when needed
            payload = mock_concept_merge_payload(vault, slug)
            self.assertTrue(payload.get("needs_review_banner"))
            # agent-like payload without banner must fail
            bad = {
                "definition": "T cell exhaustion is a dysfunctional state described across the two fixtures in this vault for testing.",
                "synthesis": [
                    {"claim": "Exhaustion appears in both papers.", "source_ids": [ZHENG_ID, REL_ID]}
                ],
                "related_concepts": [],
                "needs_review_banner": None,
            }
            with self.assertRaises(ValueError):
                apply_concept_merge(vault, slug, bad, backend="agent")
            good = dict(bad)
            good["needs_review_banner"] = "Contributing source needs_review — do not silent-merge."
            ok = apply_concept_merge(vault, slug, good, backend="agent")
            self.assertTrue(ok["ok"])
            md = (vault / "wiki" / "concepts" / f"{slug}.md").read_text(encoding="utf-8")
            self.assertIn("needs_review", md)

    def test_validate_reject_unknown_source(self):
        with self.assertRaises(ValueError):
            validate_concept_merge_payload(
                {
                    "definition": "A sufficiently long definition of the concept for schema checks here.",
                    "synthesis": [{"claim": "x", "source_ids": ["doi:unknown"]}],
                },
                allowed_source_ids={"doi:known"},
            )

    def test_redirect_merge_into(self):
        with tempfile.TemporaryDirectory() as td:
            col = Path(td) / "collections" / "c1"
            init_vault(col, collection_id="c1")
            vault = col / "wiki-vault"
            ingest_publication_to_wiki(
                vault, _zheng_pub(), compile_mode="mock", llm_json=FIX / "mock_compile_zheng.json"
            )
            ingest_publication_to_wiki(
                vault,
                _related_pub(),
                compile_mode="mock",
                llm_json=FIX / "mock_compile_related_exhaustion.json",
            )
            slug = slugify("T cell exhaustion")
            payload = mock_concept_merge_payload(vault, slug)
            payload["canonical_slug"] = "t-cell-exhaustion-canonical"
            payload["merge_into"] = "t-cell-exhaustion-canonical"
            # re-validate not needed; apply handles
            # allowed ids already in payload from mock — apply revalidates
            # strip canonical into apply via dict
            result = apply_concept_merge(vault, slug, payload, backend="mock")
            self.assertTrue(result["ok"])
            old = (vault / "wiki" / "concepts" / f"{slug}.md").read_text(encoding="utf-8")
            self.assertIn("Merged into", old)
            self.assertIn("t-cell-exhaustion-canonical", old)
            newp = vault / "wiki" / "concepts" / "t-cell-exhaustion-canonical.md"
            self.assertTrue(newp.is_file())
            self.assertIn("## Across papers", newp.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
