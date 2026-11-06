#!/usr/bin/env python3
"""Unit tests for wiki agent-compile schema + mock/file providers (v0.4)."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "lit-wiki" / "scripts"))

from wiki_llm import (  # noqa: E402
    CompileRequestNeeded,
    build_compile_request,
    compile_with_provider,
    file_compile,
    llm_payload_to_extract,
    mock_compile,
    resolve_compile_mode,
    sanitize_entity_symbol,
    validate_compile_payload,
)
from wiki_lib import compile_publication  # noqa: E402


ZHENG = {
    "id": "doi:10.1016/j.cell.2017.05.035",
    "title": "Landscape of Infiltrating T Cells in Liver Cancer Revealed by Single-Cell Sequencing",
    "abstract": (
        "Here we perform deep single-cell RNA sequencing on 5,063 single T cells from six HCC patients. "
        "We identify 11 T cell subsets. Exhausted CD8+ T cells and Tregs are enriched. "
        "LAYN (layilin) is upregulated on immunosuppressive T cells."
    ),
    "year": 2017,
    "doi": "10.1016/j.cell.2017.05.035",
    "pmid": "28622514",
}


class SchemaTests(unittest.TestCase):
    def test_validate_ok(self):
        p = validate_compile_payload(
            {
                "summary": "A grounded synthesis about HCC T cell exhaustion based on the abstract only.",
                "key_claims": ["Exhausted CD8 T cells are enriched in tumor."],
                "methods": ["scRNA-seq"],
                "concepts": [{"name": "T cell exhaustion"}],
                "entities": [{"symbol": "LAYN", "name": "layilin", "kind": "gene"}],
                "tags": ["HCC", "exhaustion"],
            }
        )
        self.assertIn("schema", p)
        ex = llm_payload_to_extract(p)
        self.assertEqual(ex["entities"], ["LAYN"])
        self.assertEqual(ex["compile_backend"], "agent")
        self.assertEqual(ex["claims"][0].startswith("Exhausted"), True)

    def test_validate_rejects_short_summary(self):
        with self.assertRaises(ValueError):
            validate_compile_payload(
                {"summary": "too short", "key_claims": ["x" * 20], "methods": [], "concepts": [], "entities": [], "tags": []}
            )

    def test_fragment_sanitize(self):
        self.assertIsNone(sanitize_entity_symbol("l1"))
        self.assertIsNone(sanitize_entity_symbol("IFN"))
        self.assertEqual(sanitize_entity_symbol("layilin"), "LAYN")
        self.assertEqual(sanitize_entity_symbol("PD-L1"), "PD-L1")


class ProviderTests(unittest.TestCase):
    def test_default_mode_is_agent(self):
        self.assertEqual(resolve_compile_mode(None), "agent")
        self.assertEqual(resolve_compile_mode("llm"), "agent")  # legacy alias
        with self.assertRaises(ValueError):
            resolve_compile_mode("openai-http")

    def test_agent_without_json_raises_request(self):
        with self.assertRaises(CompileRequestNeeded) as cm:
            compile_with_provider(ZHENG, mode="agent")
        self.assertEqual(cm.exception.request["action"], "compile_request")
        self.assertIn("summary", cm.exception.request["response_schema"]["properties"])

    def test_build_compile_request(self):
        req = build_compile_request(ZHENG, response_path="/tmp/r.json")
        self.assertIn("5,063", req["user_prompt"])
        self.assertIn("layilin", req["system_prompt"].lower() + req["user_prompt"].lower())

    def test_mock(self):
        ex = mock_compile(ZHENG)
        self.assertEqual(ex["compile_backend"], "mock")
        self.assertIn("LAYN", ex["entities"])
        self.assertNotIn("l1", [e.lower() for e in ex["entities"]])
        self.assertGreater(len(ex["summary"]), 60)
        self.assertFalse(ex["summary"].startswith("Here we perform"))

    def test_file_provider_as_agent_apply(self):
        payload = {
            "publication_id": ZHENG["id"],
            "summary": "Synthesised wiki summary of Zheng 2017 scRNA-seq atlas of HCC T cells with exhaustion and LAYN.",
            "key_claims": [
                "5,063 T cells from six HCC patients were profiled by scRNA-seq.",
                "Eleven T cell subsets were delineated with TCR information.",
                "LAYN marks immunosuppressive / activated T cell states in the abstract.",
            ],
            "methods": ["scRNA-seq", "TCR sequencing"],
            "concepts": [
                {"name": "T cell exhaustion", "slug_hint": "t-cell-exhaustion"},
                {"name": "hepatocellular carcinoma"},
            ],
            "entities": [
                {"symbol": "LAYN", "name": "layilin", "kind": "gene"},
                {"symbol": "CD8", "kind": "protein"},
            ],
            "tags": ["HCC", "exhaustion", "scRNA-seq", "LAYN"],
        }
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "compile.json"
            path.write_text(json.dumps([payload]), encoding="utf-8")
            ex = file_compile(ZHENG, path, backend="agent")
            self.assertEqual(ex["compile_backend"], "agent")
            self.assertIn("LAYN", ex["entities"])
            self.assertTrue(ex["summary"].startswith("Synthesised"))

    def test_compile_publication_mock(self):
        ex = compile_publication(ZHENG, mode="mock")
        self.assertEqual(ex["compile_backend"], "mock")
        self.assertTrue(ex["compile_label"].startswith("mock-v0.4"))


if __name__ == "__main__":
    unittest.main()
