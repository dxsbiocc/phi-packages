"""Unit and contract tests for the protein-design domain plugin:
  resources/plugins/protein-design/
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PLUGIN_DIR = REPO_ROOT / "resources" / "plugins" / "protein-design"
PREP_SCRIPTS = PLUGIN_DIR / "skills" / "structure-input-prep" / "scripts"
EVAL_SCRIPTS = PLUGIN_DIR / "skills" / "design-evaluation" / "scripts"


def _write_test_pdb(path: Path) -> None:
    lines = [
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n",
        "ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C\n",
        "ATOM      3  CA  GLY A   2       4.000   2.800   0.000  1.00 20.00           C\n",
        # Gap at residue 3-4 -> jumps to 5
        "ATOM      4  CA  SER A   5       7.600   3.700   0.500  1.00 20.00           C\n",
        "ATOM      5  CA  ASP B  10       2.800   5.600   0.000  1.00 20.00           C\n",
        "HETATM    6  C1  ATP L   1       4.000   2.800   3.200  1.00 20.00           C\n",
        "HETATM    7  O1  ATP L   1       5.000   2.800   3.200  1.00 20.00           O\n",
        "HETATM    8  O   HOH W   1      10.000  10.000  10.000  1.00 20.00           O\n",
        "END\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")


class TestProteinDesignPluginContract(unittest.TestCase):
    def test_plugin_layout_and_reserved_directories(self) -> None:
        self.assertTrue((PLUGIN_DIR / "phi-package.yaml").is_file())
        self.assertTrue((PLUGIN_DIR / "README.md").is_file())
        self.assertTrue((PLUGIN_DIR / "icon.svg").is_file())
        self.assertTrue((PLUGIN_DIR / "agents" / "ProteinDesign.md").is_file())

        # Plugin contract forbids reserved directories
        for forbidden in ("wrappers", "mcp", "orchestrator"):
            self.assertFalse((PLUGIN_DIR / forbidden).exists(), f"Reserved dir {forbidden}/ must not exist in plugin")

        manifest = (PLUGIN_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("id: protein-design", manifest)
        self.assertIn("toolPrefix: pd", manifest)
        for skill in ("protein-design-workflow", "structure-input-prep", "design-evaluation"):
            self.assertIn(f"skills/{skill}", manifest)
            self.assertTrue((PLUGIN_DIR / "skills" / skill / "SKILL.md").is_file())


class TestProteinDesignPluginScripts(unittest.TestCase):
    def test_inspect_structure_and_validate_rfd3_spec(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            pdb = tmp_dir / "target.pdb"
            _write_test_pdb(pdb)

            proc = subprocess.run(
                [sys.executable, str(PREP_SCRIPTS / "inspect_structure.py"), "--input", str(pdb)],
                capture_output=True,
                text=True,
                check=True,
            )
            info = json.loads(proc.stdout)
            self.assertEqual(info["num_chains"], 2)
            chain_a = next(c for c in info["chains"] if c["chain_id"] == "A")
            self.assertEqual(chain_a["gaps"], ["3-4"])
            self.assertEqual(chain_a["sequence"], "AGS")
            self.assertEqual(len(info["ligands"]), 1)
            self.assertEqual(info["ligands"][0]["resname"], "ATP")

            # Valid contig + ligand
            ok_proc = subprocess.run(
                [
                    sys.executable,
                    str(PREP_SCRIPTS / "validate_rfd3_spec.py"),
                    "--contig",
                    "A1-5/0 60-80",
                    "--input-structure",
                    str(pdb),
                    "--ligand",
                    "ATP",
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            ok_res = json.loads(ok_proc.stdout)
            self.assertTrue(ok_res["valid"], ok_res)

            # Invalid chain + missing ligand
            bad_proc = subprocess.run(
                [
                    sys.executable,
                    str(PREP_SCRIPTS / "validate_rfd3_spec.py"),
                    "--contig",
                    "Z1-10/0 60-80",
                    "--input-structure",
                    str(pdb),
                    "--ligand",
                    "HEM",
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            bad_res = json.loads(bad_proc.stdout)
            self.assertFalse(bad_res["valid"])
            self.assertGreaterEqual(len(bad_res["errors"]), 2)

    def test_funnel_report(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            merged = tmp_dir / "merged.tsv"
            filtered = tmp_dir / "filtered.tsv"
            summary = tmp_dir / "summary.json"

            merged.write_text(
                "design_id\tparent_design_id\tcomplex_plddt\tiptm\tca_rmsd\tpassed_filter\n"
                "d1\tp1\t0.88\t0.81\t1.2\ttrue\n"
                "d2\tp1\t0.84\t0.76\t1.5\ttrue\n"
                "d3\tp2\t0.91\t0.85\t1.0\ttrue\n"
                "d4\tp3\t0.62\t0.41\t3.4\tfalse\n",
                encoding="utf-8",
            )
            filtered.write_text(
                "design_id\tparent_design_id\tcomplex_plddt\tiptm\tca_rmsd\tpassed_filter\n"
                "d3\tp2\t0.91\t0.85\t1.0\ttrue\n"
                "d1\tp1\t0.88\t0.81\t1.2\ttrue\n",
                encoding="utf-8",
            )
            summary.write_text(
                json.dumps({"expr": "complex_plddt >= 0.8 and iptm >= 0.7", "rank_by": "iptm"}),
                encoding="utf-8",
            )

            proc = subprocess.run(
                [
                    sys.executable,
                    str(EVAL_SCRIPTS / "funnel_report.py"),
                    "--merged-tsv",
                    str(merged),
                    "--filtered-tsv",
                    str(filtered),
                    "--summary-json",
                    str(summary),
                    "--top-k",
                    "2",
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            rep = json.loads(proc.stdout)
            self.assertEqual(rep["total_designs"], 4)
            self.assertEqual(rep["passed_filter"], 3)
            self.assertEqual(rep["selected_designs"], 2)
            self.assertAlmostEqual(rep["pass_rate"], 0.75, places=3)
            self.assertEqual(rep["unique_parents_total"], 3)
            self.assertEqual(rep["unique_parents_selected"], 2)
            self.assertEqual(len(rep["top_candidates"]), 2)
            self.assertEqual(rep["top_candidates"][0]["design_id"], "d3")
            self.assertIn("iptm", rep["metric_quantiles_all"])


if __name__ == "__main__":
    unittest.main()
