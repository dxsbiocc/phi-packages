#!/usr/bin/env python3
"""End-to-end checks of the lit_* tool scripts, run as the engine runs them (offline)."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "skills" / "manage-library" / "scripts"
FX = ROOT / "tests" / "fixtures"


def run(script: str, *args: str, ok: bool = True) -> dict:
    proc = subprocess.run(
        [sys.executable, str(SCRIPTS / script), *args], capture_output=True, text=True, check=False
    )
    lines = proc.stdout.strip().splitlines()
    assert lines, f"{script} printed nothing: {proc.stderr[-400:]}"
    out = json.loads(lines[-1])
    assert out.get("ok") is ok, f"{script} {args}: {out}\n{proc.stderr[-400:]}"
    return out


class Project(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.addCleanup(self._td.cleanup)
        self.root = Path(self._td.name)
        self.p = ["--project_root", str(self.root)]

    def make_collection(self, name="Demo") -> str:
        return run("collection.py", "--action", "create", *self.p, "--name", name)["collection_id"]

    def search(self, query="exhaustion", *extra):
        return run("search.py", *self.p, "--query", query, "--mock", *extra)


class StatusAndSearchTests(Project):
    def test_status_before_and_after(self):
        out = run("status.py", *self.p)
        self.assertFalse(out["initialized"])
        self.make_collection()
        out = run("status.py", *self.p)
        self.assertEqual(out["collections"][0]["n_papers"], 0)
        self.assertTrue(any("lit_search" in s for s in out["collections"][0]["next_steps"]))

    def test_search_uses_the_only_collection_and_records_the_run(self):
        cid = self.make_collection()
        out = self.search("t cell exhaustion", "--sources", "pubmed,europepmc")
        self.assertEqual(out["collection_id"], cid)
        self.assertEqual([s["source"] for s in out["per_source"]], ["pubmed", "europepmc"])
        self.assertGreaterEqual(out["ingest"]["counts"]["inserted"], 1)
        run_json = json.loads((Path(out["hits_path"]).parent / "run.json").read_text())
        self.assertEqual(run_json["query"], "t cell exhaustion")
        self.assertNotIn("gold_set", run_json)

    def test_search_rejects_unsupported_sources_and_too_many(self):
        self.make_collection()
        bad = run("search.py", *self.p, "--query", "x", "--sources", "pmc", "--mock", ok=False)
        self.assertIn("unsupported", bad["error"])
        many = run("search.py", *self.p, "--query", "x", "--mock",
                   "--sources", "pubmed,europepmc,openalex,crossref,arxiv", ok=False)
        self.assertIn("at most", many["error"])

    def test_several_collections_require_an_id(self):
        self.make_collection("A")
        self.make_collection("B")
        out = run("search.py", *self.p, "--query", "x", "--mock", ok=False)
        self.assertIn("several collections", out["error"])


class ImportQueryGetExportTests(Project):
    def setUp(self):
        super().setUp()
        self.cid = self.make_collection()
        self.search()
        run("ingest.py", *self.p, "--input_path", str(FX / "sample.bib"), "--no_enrich")

    def test_query_finds_papers_and_flags_wiki_pages(self):
        out = run("query.py", *self.p, "--query", "mock")
        self.assertGreaterEqual(out["n"], 1)
        self.assertFalse(out["papers"][0]["has_wiki_page"])

    def test_get_by_doi_and_unknown(self):
        out = run("get.py", *self.p, "--publication_id", "10.1000/mock.lit.001")
        self.assertEqual(out["publication"]["id"], "doi:10.1000/mock.lit.001")
        self.assertEqual(out["collections"][0]["id"], self.cid)
        bad = run("get.py", *self.p, "--publication_id", "10.9999/nope", ok=False)
        self.assertIn("no saved paper", bad["error"])

    def test_export_formats(self):
        for fmt, marker in (("bibtex", "@article{"), ("ris", "TY  - JOUR"), ("csv", "cite_key,"), ("search_log", '"publication_id"')):
            out = run("export.py", *self.p, "--format", fmt)
            self.assertIn(marker, Path(out["out_path"]).read_text(encoding="utf-8"), fmt)
        csl = run("export.py", *self.p, "--format", "csl-json")
        items = json.loads(Path(csl["out_path"]).read_text(encoding="utf-8"))
        self.assertTrue(all("id" in i and "title" in i for i in items))
        self.assertEqual(len({i["id"] for i in items}), len(items))

    def test_import_folder_records_files_and_access_origin(self):
        folder = self.root / "papers"
        folder.mkdir()
        shutil.copy(FX / "tiny.pdf", folder / "tiny.pdf")
        shutil.copy(FX / "sample.ris", folder / "sample.ris")
        out = run("ingest.py", *self.p, "--input_path", str(folder), "--no_enrich", "--access_origin", "institutional")
        self.assertGreaterEqual(out["n_parsed"], 2)
        got = run("get.py", *self.p, "--publication_id", "10.1000/tiny.pdf.001")
        self.assertEqual(got["files"][0]["access_origin"], "institutional")


class CurateTests(Project):
    def test_queue_and_merge_dry_run(self):
        self.make_collection()
        self.search()
        self.assertTrue(run("curate.py", "--action", "queue", *self.p)["backlog_empty"])
        self.assertEqual(run("curate.py", "--action", "merge_ids", *self.p)["mode"], "dry-run")
        run("curate.py", "--action", "rebind", *self.p, ok=False)  # needs scan_dir


class WikiBatchTests(Project):
    def setUp(self):
        super().setUp()
        self.cid = self.make_collection()
        self.search()

    def payload(self, pid, **over):
        base = {
            "publication_id": pid,
            "summary": "This paper reports a finding about T cell exhaustion, grounded in the abstract.",
            "key_claims": ["A claim from the abstract"],
            "methods": ["flow cytometry"],
            "concepts": [{"name": "T cell exhaustion", "evidence": "named"}, {"name": "PD-1"}],
            "tags": ["exhaustion"],
        }
        base.update(over)
        return base

    def test_batch_roundtrip_with_a_bad_result(self):
        pend = run("wiki_pending.py", *self.p)
        self.assertEqual(pend["n_in_batch"], 2)
        req = json.loads(Path(pend["request_path"]).read_text())
        ids = [i["publication_id"] for i in req["items"]]
        results = [self.payload(i) for i in ids] + [self.payload("doi:unknown"), self.payload(ids[0], summary="short")]
        Path(pend["response_path"]).write_text(json.dumps({"results": results}))
        out = run("wiki_apply.py", *self.p)
        self.assertEqual(out["n_written"], 2)
        self.assertEqual(out["n_failed"], 2)
        self.assertTrue(all(w["tags_filled"] for w in out["written"]))
        # two papers mention the same concepts, so both need a Definition
        self.assertEqual({c["slug"] for c in out["concepts_needing_definition"]}, {"t-cell-exhaustion", "pd-1"})
        self.assertEqual(run("wiki_pending.py", *self.p)["n_in_batch"], 0)
        lint = run("wiki_lint.py", *self.p)
        self.assertEqual(lint["n_errors"], 0)
        self.assertEqual(lint["n_info"], 2)

    def test_defining_a_concept_clears_the_lint_note(self):
        pend = run("wiki_pending.py", *self.p)
        ids = [i["publication_id"] for i in json.loads(Path(pend["request_path"]).read_text())["items"]]
        Path(pend["response_path"]).write_text(json.dumps({"results": [self.payload(i) for i in ids]}))
        out = run("wiki_apply.py", *self.p)
        page = Path(pend["response_path"]).parents[1] / "wiki" / "concepts" / "pd-1.md"
        text = page.read_text()
        page.write_text(text.replace("_Not written yet. Write it once two or more papers mention this concept._", "A checkpoint receptor."))
        issues = run("wiki_lint.py", *self.p)["issues"]
        slugs = {Path(i["path"]).stem for i in issues if i["code"] == "concept_needs_definition"}
        self.assertEqual(slugs, {"t-cell-exhaustion"})
        self.assertEqual(out["n_failed"], 0)

    def test_query_marks_papers_with_pages(self):
        pend = run("wiki_pending.py", *self.p)
        ids = [i["publication_id"] for i in json.loads(Path(pend["request_path"]).read_text())["items"]]
        Path(pend["response_path"]).write_text(json.dumps({"results": [self.payload(i) for i in ids]}))
        run("wiki_apply.py", *self.p)
        out = run("query.py", *self.p, "--query", "exhaustion")
        self.assertTrue(all(p["has_wiki_page"] for p in out["papers"]))
        self.assertTrue(out["wiki_hits"])

    def test_lint_flags_a_paper_missing_from_the_index(self):
        pend = run("wiki_pending.py", *self.p)
        ids = [i["publication_id"] for i in json.loads(Path(pend["request_path"]).read_text())["items"]]
        Path(pend["response_path"]).write_text(json.dumps({"results": [self.payload(ids[0])]}))
        run("wiki_apply.py", *self.p)
        index = Path(pend["response_path"]).parents[1] / "wiki" / "index.md"
        index.write_text("# Index\n")
        codes = {i["code"] for i in run("wiki_lint.py", *self.p)["issues"]}
        self.assertIn("missing_index_entry", codes)


if __name__ == "__main__":
    unittest.main()
