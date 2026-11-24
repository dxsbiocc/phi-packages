"""Unit, contract, parameter parity, and relative include resolution tests for
AlphaFold 3 and Phase 3 protein design subworkflows and end-to-end workflow:
  - resources/wrappers/modules/local/alphafold3
  - resources/wrappers/subworkflows/local/design_postprocess
  - resources/wrappers/subworkflows/local/seq_fold_eval
  - resources/wrappers/workflows/local/protein_binder_design
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
WRAPPERS_BASE = REPO_ROOT / "resources" / "wrappers"
AF3_DIR = WRAPPERS_BASE / "modules" / "local" / "alphafold3"


class TestAlphaFold3Scripts(unittest.TestCase):
    def test_prepare_and_compile_alphafold3_from_sequences_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            seq_tsv = tmp_dir / "sequences.tsv"
            seq_tsv.write_text(
                "design_id\tparent_design_id\tsequence\tparent_structure_path\n"
                "BB_0_mpnn_1\tBB_0\tMKTAYIAKQRQ\t/tmp/bb_0.pdb\n",
                encoding="utf-8",
            )

            af3_in = tmp_dir / "af3_inputs"
            af3_manifest = tmp_dir / "af3_manifest.json"
            subprocess.run(
                [
                    sys.executable,
                    str(AF3_DIR / "bin" / "prepare_af3_inputs.py"),
                    "--sequences-tsv",
                    str(seq_tsv),
                    "--target-sequence",
                    "ACDEFGHIKLMNPQRSTVWY",
                    "--ligand-ccd",
                    "ATP",
                    "--model-seeds",
                    "1,2",
                    "--out-dir",
                    str(af3_in),
                    "--out-manifest",
                    str(af3_manifest),
                ],
                check=True,
            )

            job_json_path = af3_in / "bb_0_mpnn_1.json"
            self.assertTrue(job_json_path.is_file())
            job_data = json.loads(job_json_path.read_text(encoding="utf-8"))
            self.assertEqual(job_data["dialect"], "alphafold3")
            self.assertEqual(job_data["version"], 1)
            self.assertEqual(job_data["modelSeeds"], [1, 2])
            self.assertEqual(len(job_data["sequences"]), 3)
            self.assertEqual(job_data["sequences"][0]["protein"]["unpairedMsa"], "")
            self.assertEqual(job_data["sequences"][2]["ligand"]["ccdCodes"], ["ATP"])

            # Simulate AF3 raw output directory
            job_out = tmp_dir / "af3_raw" / "bb_0_mpnn_1"
            job_out.mkdir(parents=True)
            (job_out / "bb_0_mpnn_1_model.cif").write_text("data_af3\n#\n", encoding="utf-8")
            (job_out / "bb_0_mpnn_1_summary_confidences.json").write_text(
                json.dumps(
                    {
                        "ranking_score": 0.87,
                        "ptm": 0.85,
                        "iptm": 0.82,
                        "fraction_disordered": 0.04,
                        "has_clash": 0.0,
                        "chain_pair_pae_min": [[0.4, 3.1], [2.4, 0.5]],
                    }
                ),
                encoding="utf-8",
            )
            (job_out / "bb_0_mpnn_1_confidences.json").write_text(
                json.dumps({"atom_plddts": [88.0, 92.0, 90.0]}),
                encoding="utf-8",
            )

            flat_structs = tmp_dir / "af3_structs"
            out_tsv = tmp_dir / "af3.tsv"
            out_json = tmp_dir / "af3.json"
            subprocess.run(
                [
                    sys.executable,
                    str(AF3_DIR / "bin" / "compile_af3_outputs.py"),
                    "--out-dir",
                    str(tmp_dir / "af3_raw"),
                    "--manifest",
                    str(af3_manifest),
                    "--out-structures-dir",
                    str(flat_structs),
                    "--out-tsv",
                    str(out_tsv),
                    "--out-json",
                    str(out_json),
                ],
                check=True,
            )

            rows = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["design_id"], "BB_0_mpnn_1")
            self.assertEqual(row["parent_design_id"], "BB_0")
            self.assertAlmostEqual(row["confidence_score"], 0.87, places=3)
            self.assertAlmostEqual(row["ranking_score"], 0.87, places=3)
            self.assertAlmostEqual(row["plddt"], 90.0, places=3)
            self.assertAlmostEqual(row["complex_plddt"], 0.90, places=3)
            self.assertAlmostEqual(row["min_pae_interaction"], 2.4, places=3)
            self.assertTrue((flat_structs / "BB_0_mpnn_1.cif").is_file())


class TestProteinDesignComposites(unittest.TestCase):
    def _assert_includes_resolve(self, nf_path: Path) -> None:
        text = nf_path.read_text(encoding="utf-8")
        matches = list(re.finditer(r"include\s*\{[^}]+\}\s*from\s*['\"]([^'\"]+)['\"]", text))
        self.assertGreater(len(matches), 0, f"Expected at least one include in {nf_path}")
        for match in matches:
            rel = match.group(1)
            target = (nf_path.parent / rel).resolve()
            if not target.is_file():
                target_main = (nf_path.parent / rel / "main.nf").resolve()
                self.assertTrue(
                    target_main.is_file(),
                    f"Broken include in {nf_path}: {rel} (neither {target} nor {target_main} exists)",
                )

    def _assert_yaml_and_params_parity(self, wrapper_dir: Path, expected_id: str) -> None:
        yaml_text = (wrapper_dir / "wrapper.yaml").read_text(encoding="utf-8")
        params_json = json.loads((wrapper_dir / "params.json").read_text(encoding="utf-8"))
        self.assertTrue((wrapper_dir / "main.nf").is_file())
        self.assertTrue((wrapper_dir / "nextflow.config").is_file())
        self._assert_includes_resolve(wrapper_dir / "main.nf")

        self.assertIn(f"id: {expected_id}", yaml_text)
        self.assertIn("primary: true", yaml_text)
        self.assertNotIn("\n    default:", yaml_text, "defaults must live only in wrapper/params.json")
        self.assertIn("outdir", params_json)
        self.assertIsInstance(params_json, dict)

        in_params = False
        yaml_param_keys: set[str] = set()
        for line in yaml_text.splitlines():
            if line.startswith("params:"):
                in_params = True
                continue
            if in_params and line and not line.startswith(" "):
                break
            if in_params:
                m = re.match(r"^  ([a-z][a-z0-9_]*):\s*$", line)
                if m:
                    yaml_param_keys.add(m.group(1))

        self.assertEqual(
            yaml_param_keys,
            set(params_json.keys()),
            f"Param mismatch between wrapper.yaml and params.json in {wrapper_dir}",
        )

    def test_all_six_modules_param_parity_and_includes(self) -> None:
        modules = [
            ("protein_metrics", "local/modules/protein-metrics"),
            ("design_filter", "local/modules/design-filter"),
            ("rfd3", "local/modules/rfd3"),
            ("ligandmpnn", "local/modules/ligandmpnn"),
            ("boltz2", "local/modules/boltz2"),
            ("alphafold3", "local/modules/alphafold3"),
        ]
        for folder, expected_id in modules:
            self._assert_yaml_and_params_parity(
                WRAPPERS_BASE / "modules" / "local" / folder / "wrapper",
                expected_id,
            )

    def test_subworkflows_and_workflow_contract_parity_and_includes(self) -> None:
        composites = [
            (WRAPPERS_BASE / "subworkflows" / "local" / "design_postprocess", "local/subworkflows/design-postprocess", True),
            (WRAPPERS_BASE / "subworkflows" / "local" / "seq_fold_eval", "local/subworkflows/seq-fold-eval", True),
            (WRAPPERS_BASE / "workflows" / "local" / "protein_binder_design", "local/workflows/protein-binder-design", False),
        ]
        for comp_dir, expected_id, require_meta in composites:
            main_nf = comp_dir / "main.nf"
            self.assertTrue(main_nf.is_file(), f"missing main.nf in {comp_dir}")
            if require_meta:
                self.assertTrue((comp_dir / "meta.yml").is_file(), f"missing meta.yml in {comp_dir}")
            dag_path = comp_dir / "wrapper" / "dag.mmd"
            self.assertTrue(dag_path.is_file(), f"missing wrapper/dag.mmd in {comp_dir}")
            self.assertIn("flowchart", dag_path.read_text(encoding="utf-8"))
            self._assert_includes_resolve(main_nf)
            self._assert_yaml_and_params_parity(comp_dir / "wrapper", expected_id)

        # Verify seq_fold_eval wires both BOLTZ2 and ALPHAFOLD3
        sfe_text = (WRAPPERS_BASE / "subworkflows" / "local" / "seq_fold_eval" / "main.nf").read_text(encoding="utf-8")
        self.assertIn("BOLTZ2", sfe_text)
        self.assertIn("ALPHAFOLD3", sfe_text)

    def test_design_postprocess_foldseek_cluster_wiring_and_bridge(self) -> None:
        dp_dir = WRAPPERS_BASE / "subworkflows" / "local" / "design_postprocess"
        dp_nf = (dp_dir / "main.nf").read_text(encoding="utf-8")
        dp_params = json.loads((dp_dir / "wrapper" / "params.json").read_text(encoding="utf-8"))

        # Verify subworkflow DSL wiring and auto-defaulting of group_by=cluster_id
        self.assertIn("include { FOLDSEEK", dp_nf)
        self.assertIn("run_foldseek_cluster", dp_nf)
        self.assertIn("run_foldseek_cluster ? 'cluster_id' : ''", dp_nf)
        self.assertIn("foldseek_tsv = ch_foldseek_tsv", dp_nf)
        self.assertFalse(dp_params["run_foldseek_cluster"])
        self.assertIn("foldseek_coverage", dp_params)
        self.assertIn("foldseek_tmscore_threshold", dp_params)

        # End-to-end Python bridge: Foldseek easy-cluster -> cluster_id -> DESIGN_FILTER diversity cap
        foldseek_script = WRAPPERS_BASE / "modules" / "local" / "foldseek" / "bin" / "compile_foldseek_outputs.py"
        filter_script = WRAPPERS_BASE / "modules" / "local" / "design_filter" / "bin" / "merge_and_filter_designs.py"

        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            raw_clu = tmp_dir / "designs_clu.tsv"
            raw_clu.write_text(
                "D1.cif\tD1.cif\n"
                "D1.cif\tD2.cif\n"
                "D1.cif\tD3.cif\n"
                "D4.cif\tD4.cif\n",
                encoding="utf-8",
            )
            foldseek_tsv = tmp_dir / "foldseek_cluster.tsv"
            foldseek_json = tmp_dir / "foldseek_cluster.json"
            subprocess.run(
                [
                    sys.executable,
                    str(foldseek_script),
                    "--mode",
                    "cluster",
                    "--input-tsv",
                    str(raw_clu),
                    "--out-tsv",
                    str(foldseek_tsv),
                    "--out-json",
                    str(foldseek_json),
                ],
                check=True,
            )

            metrics_tsv = tmp_dir / "metrics.tsv"
            metrics_tsv.write_text(
                "design_id\tconfidence_score\tclash_count\n"
                "D1\t0.72\t0\n"
                "D2\t0.91\t0\n"
                "D3\t0.85\t2\n"
                "D4\t0.80\t0\n",
                encoding="utf-8",
            )

            merged_tsv = tmp_dir / "merged.tsv"
            filtered_tsv = tmp_dir / "filtered.tsv"
            summary_json = tmp_dir / "summary.json"
            subprocess.run(
                [
                    sys.executable,
                    str(filter_script),
                    "--tables",
                    str(metrics_tsv),
                    str(foldseek_tsv),
                    "--id-column",
                    "design_id",
                    "--join-mode",
                    "outer",
                    "--expr",
                    "clash_count == 0 and confidence_score >= 0.7",
                    "--rank-by",
                    "confidence_score",
                    "--rank-direction",
                    "max",
                    "--group-by",
                    "cluster_id",
                    "--max-per-group",
                    "1",
                    "--out-merged",
                    str(merged_tsv),
                    "--out-filtered",
                    str(filtered_tsv),
                    "--out-summary",
                    str(summary_json),
                ],
                check=True,
            )

            summary = json.loads(summary_json.read_text(encoding="utf-8"))
            self.assertEqual(summary["total_designs"], 4)
            self.assertEqual(summary["passed_filter"], 3)
            self.assertEqual(summary["selected_designs"], 2)
            self.assertEqual(summary["group_by"], "cluster_id")
            self.assertEqual(summary["max_per_group"], 1)

            filtered_lines = [
                line.split("\t")
                for line in filtered_tsv.read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            header = filtered_lines[0]
            id_idx = header.index("design_id")
            clu_idx = header.index("cluster_id")
            selected_ids = [row[id_idx] for row in filtered_lines[1:]]
            selected_clusters = [row[clu_idx] for row in filtered_lines[1:]]
            # D2 (0.91) beats D1 (0.72) inside clu_1; D4 (0.80) wins clu_2
            self.assertEqual(selected_ids, ["D2", "D4"])
            self.assertEqual(selected_clusters, ["clu_1", "clu_2"])


if __name__ == "__main__":
    unittest.main()
