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
sys.path.insert(0, str(ROOT / "skills" / "lit-wiki" / "scripts"))
sys.path.insert(0, str(ROOT / "skills" / "literature-library" / "scripts"))
sys.path.insert(0, str(ROOT / "skills" / "lit-wiki" / "scripts"))

from wiki_lib import (  # noqa: E402
    init_vault,
    lint_vault,
    resolve_related_concept_slug,
    resolve_related_concepts,
    slugify,
)
from wiki_llm import apply_concept_merge, validate_compile_payload  # noqa: E402
from db import (  # noqa: E402
    coverage_actionable_hints,
    init_db,
    open_db,
    record_pub_file,
    utc_now,
)
from rebind_files import rebind_pub_files  # noqa: E402
from wiki_agent_status import agent_status  # noqa: E402


FIX = ROOT / "tests" / "fixtures"


class RelatedConceptResolveTests(unittest.TestCase):
    def test_nsclc_alias_no_broken_link(self):
        with tempfile.TemporaryDirectory() as td:
            col = Path(td) / "collections" / "c1"
            init_vault(col, collection_id="c1", collection_name="ux")
            vault = col / "wiki-vault"
            concepts = vault / "wiki" / "concepts"
            concepts.mkdir(parents=True, exist_ok=True)
            # Existing long-name page
            (concepts / "non-small-cell-lung-cancer.md").write_text(
                "# Non-small cell lung cancer\n\nSlug: `non-small-cell-lung-cancer`\n\n"
                "## Definition\n\nPlaceholder definition for alias resolve tests.\n\n"
                "## Across papers\n\n- stub (`doi:x`)\n\n"
                "## Sources\n\n- [stub](../sources/x.md)\n\n"
                "## Related concepts\n\n",
                encoding="utf-8",
            )
            (concepts / "t-cell-exhaustion.md").write_text(
                "# T cell exhaustion\n\nSlug: `t-cell-exhaustion`\n\n"
                "## Definition\n\nPlaceholder.\n\n## Across papers\n\n- a (`doi:a`)\n\n"
                "## Sources\n\n- [a](../sources/a.md)\n\n## Related concepts\n\n",
                encoding="utf-8",
            )
            # contribs so apply_concept_merge validation can cite
            contrib = concepts / "_contrib" / "t-cell-exhaustion"
            contrib.mkdir(parents=True)
            for sid, title in (("doi_a", "Paper A"), ("doi_b", "Paper B")):
                (contrib / f"{sid}.json").write_text(
                    json.dumps(
                        {
                            "concept_name": "T cell exhaustion",
                            "slug": "t-cell-exhaustion",
                            "source_id": f"doi:{sid}",
                            "source_safe_id": sid,
                            "source_title": title,
                            "source_rel": f"../sources/{sid}.md",
                            "summary_snippet": "Exhausted T cells in tumors.",
                            "claims": ["Exhaustion markers rise."],
                            "evidence": "exhausted T cells",
                            "needs_review": False,
                        }
                    ),
                    encoding="utf-8",
                )

            self.assertEqual(
                resolve_related_concept_slug(vault, "nsclc"),
                "non-small-cell-lung-cancer",
            )
            resolved, unresolved = resolve_related_concepts(
                vault, ["nsclc", "totally-missing-concept-xyz", "T cell exhaustion"]
            )
            slugs = {r["slug"] for r in resolved}
            self.assertIn("non-small-cell-lung-cancer", slugs)
            self.assertIn("t-cell-exhaustion", slugs)
            self.assertIn("totally-missing-concept-xyz", unresolved)

            payload = {
                "definition": (
                    "T cell exhaustion is a dysfunctional state discussed across "
                    "the two contribution fixtures in this offline test."
                ),
                "synthesis": [
                    {"claim": "Markers rise in tumors.", "source_ids": ["doi:doi_a"]},
                    {"claim": "Shared program across papers.", "source_ids": ["doi:doi_b"]},
                ],
                "related_concepts": ["nsclc", "ghost-concept-absent"],
                "aliases": [],
            }
            result = apply_concept_merge(vault, "t-cell-exhaustion", payload, backend="agent")
            self.assertTrue(result["ok"])
            self.assertTrue(any(r["slug"] == "non-small-cell-lung-cancer" for r in result["related_resolved"]))
            self.assertIn("ghost-concept-absent", result["related_unresolved"])
            page = (concepts / "t-cell-exhaustion.md").read_text(encoding="utf-8")
            self.assertIn("non-small-cell-lung-cancer.md", page)
            self.assertNotIn("](nsclc.md)", page)
            self.assertNotIn("ghost-concept-absent.md", page)
            lint = lint_vault(vault)
            broken = [i for i in lint.get("issues") or [] if i.get("code") == "broken_link"]
            # related link must not introduce broken nsclc
            self.assertFalse(any("nsclc.md" in str(i.get("link")) for i in broken), broken)


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

            q_script = ROOT / "skills" / "literature-library" / "scripts" / "query.py"
            wq_script = ROOT / "skills" / "lit-wiki" / "scripts" / "wiki_query.py"
            for script, extra in (
                (q_script, []),
                (wq_script, ["--collection_id", cid]),
            ):
                cmd = [
                    sys.executable,
                    str(script),
                    "--project_root",
                    str(root),
                    "--query",
                    "exhaustion",
                    *extra,
                ]
                if script == q_script:
                    cmd.extend(["--collection_id", cid])
                proc = subprocess.run(cmd, capture_output=True, text=True)
                self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
                out = json.loads(proc.stdout)
                self.assertTrue(out.get("ok"))


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


class CoverageHintsTests(unittest.TestCase):
    def test_hints_when_recall_low(self):
        hints = coverage_actionable_hints(
            {
                "recall_evidence": {
                    "gold_set_recall": {
                        "gold_set_n": 7,
                        "n_found": 2,
                        "n_missed": 5,
                        "recall": 2 / 7,
                        "evidence_weak": False,
                    },
                    "out_of_window": {"n_controls": 0},
                    "declared_but_uncontacted_sources": {"n_uncontacted": 0},
                },
                "collection": {"n_unenriched_shells": 0, "n_needs_review": 0},
                "hits": {"n": 8},
            }
        )
        self.assertTrue(any("gold_set recall" in h for h in hints))
        self.assertTrue(any("seed" in h.lower() or "ids" in h.lower() for h in hints))
        empty = coverage_actionable_hints(
            {
                "recall_evidence": {
                    "gold_set_recall": {"gold_set_n": 0, "recall": None},
                    "out_of_window": {},
                    "declared_but_uncontacted_sources": {},
                },
                "collection": {},
                "hits": {},
            }
        )
        self.assertEqual(empty, [])


class AgentStatusTests(unittest.TestCase):
    def test_agent_status_basic(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            db = root / "library.sqlite"
            init_db(db)
            cid = "col-status"
            col = root / "collections" / cid
            col.mkdir(parents=True)
            conn = open_db(db)
            conn.execute(
                "INSERT INTO collections(id,name,description,created_at,path) VALUES (?,?,?,?,?)",
                (cid, "S", "", utc_now(), str(col)),
            )
            conn.execute(
                "INSERT INTO publications(id,title,abstract,created_at,updated_at) VALUES (?,?,?,?,?)",
                (
                    "doi:10.status/1",
                    "Landscape of infiltrating T cells",
                    "Deep single-cell RNA sequencing of T cells.",
                    utc_now(),
                    utc_now(),
                ),
            )
            conn.execute(
                "INSERT INTO collection_pubs(collection_id,publication_id,added_at) VALUES (?,?,?)",
                (cid, "doi:10.status/1", utc_now()),
            )
            conn.commit()
            conn.close()
            init_vault(col, collection_id=cid, collection_name="S")
            out = agent_status(project_root=root, collection_id=cid, limit=10)
            self.assertTrue(out["ok"])
            self.assertEqual(out["n_missing_wiki_sources"], 1)
            self.assertEqual(out["missing_wiki_sources"][0]["publication_id"], "doi:10.status/1")
            self.assertTrue(out["suggested_next"])


class SoftEntityValidateTests(unittest.TestCase):
    def test_strip_implied_entities(self):
        payload = validate_compile_payload(
            {
                "summary": "A sufficiently long summary for validation of compile payload schema.",
                "key_claims": ["Claim one is grounded in the abstract text provided."],
                "methods": ["scRNA-seq"],
                "concepts": [{"name": "T cell exhaustion"}],
                "entities": [
                    {"symbol": "LAYN", "name": "layilin", "kind": "gene"},
                    {
                        "symbol": "FOXP3",
                        "name": "implied regulator",
                        "note": "not claimed in abstract",
                    },
                ],
                "tags": ["exhaustion"],
            }
        )
        syms = [e["symbol"] for e in payload["entities"]]
        self.assertIn("LAYN", syms)
        self.assertNotIn("FOXP3", syms)
        self.assertTrue(payload.get("entity_warnings"))


if __name__ == "__main__":
    unittest.main()
