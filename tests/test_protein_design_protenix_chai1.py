"""Unit and contract tests for Protenix and Chai-1 Nextflow modules and their
integration into seq_fold_eval and protein_binder_design:
  - resources/wrappers/modules/local/protenix
  - resources/wrappers/modules/local/chai1
  - resources/wrappers/subworkflows/local/seq_fold_eval
  - resources/wrappers/workflows/local/protein_binder_design
"""

from __future__ import annotations

import json
import re
import struct
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
WRAPPERS_BASE = REPO_ROOT / "resources" / "wrappers"
PROTENIX_DIR = WRAPPERS_BASE / "modules" / "local" / "protenix"
CHAI1_DIR = WRAPPERS_BASE / "modules" / "local" / "chai1"


def _make_npy_f32(value: float) -> bytes:
    header_str = "{'descr': '<f4', 'fortran_order': False, 'shape': (1,), }"
    # Pad header to 16-byte alignment (including 10-byte prefix + trailing newline)
    pad_len = (16 - ((10 + len(header_str) + 1) % 16)) % 16
    header_bytes = (header_str + (" " * pad_len) + "\n").encode("latin1")
    return b"\x93NUMPY\x01\x00" + struct.pack("<H", len(header_bytes)) + header_bytes + struct.pack("<f", value)


def _make_npy_bool(value: bool) -> bytes:
    header_str = "{'descr': '|b1', 'fortran_order': False, 'shape': (1,), }"
    pad_len = (16 - ((10 + len(header_str) + 1) % 16)) % 16
    header_bytes = (header_str + (" " * pad_len) + "\n").encode("latin1")
    return b"\x93NUMPY\x01\x00" + struct.pack("<H", len(header_bytes)) + header_bytes + bytes([1 if value else 0])


def _write_chai_npz(path: Path, agg: float, ptm: float, iptm: float, clash: bool) -> None:
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("aggregate_score.npy", _make_npy_f32(agg))
        zf.writestr("ptm.npy", _make_npy_f32(ptm))
        zf.writestr("iptm.npy", _make_npy_f32(iptm))
        zf.writestr("has_inter_chain_clashes.npy", _make_npy_bool(clash))


class TestProtenixScripts(unittest.TestCase):
    def test_prepare_and_compile_protenix(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            seq_tsv = tmp_dir / "sequences.tsv"
            seq_tsv.write_text(
                "design_id\tparent_design_id\tsequence\tparent_structure_path\n"
                "bb1_mpnn_1\tbb1\tMKTAYIAKQRQ\t/tmp/bb1.pdb\n",
                encoding="utf-8",
            )
            staged_json = tmp_dir / "protenix.json"
            manifest = tmp_dir / "manifest.json"
            subprocess.run(
                [
                    sys.executable,
                    str(PROTENIX_DIR / "bin" / "prepare_protenix_inputs.py"),
                    "--sequences-tsv",
                    str(seq_tsv),
                    "--target-sequence",
                    "ACDEFGHIKLMNPQRSTVWY",
                    "--ligand-ccd",
                    "ATP",
                    "--out-json",
                    str(staged_json),
                    "--out-manifest",
                    str(manifest),
                ],
                check=True,
            )
            jobs = json.loads(staged_json.read_text(encoding="utf-8"))
            self.assertEqual(len(jobs), 1)
            self.assertEqual(jobs[0]["name"], "bb1_mpnn_1")
            self.assertEqual(jobs[0]["sequences"][2]["ligand"]["ligand"], "CCD_ATP")

            pred_dir = tmp_dir / "protenix_raw" / "bb1_mpnn_1" / "seed_101" / "predictions"
            pred_dir.mkdir(parents=True)
            (pred_dir / "bb1_mpnn_1_sample_0.cif").write_text("data_s0\n#\n", encoding="utf-8")
            (pred_dir / "bb1_mpnn_1_sample_1.cif").write_text("data_s1\n#\n", encoding="utf-8")
            (pred_dir / "bb1_mpnn_1_summary_confidence_sample_0.json").write_text(
                json.dumps({"ranking_score": 0.71, "ptm": 0.72, "iptm": 0.68, "plddt": 81.0, "gpde": 1.8, "has_clash": False}),
                encoding="utf-8",
            )
            (pred_dir / "bb1_mpnn_1_summary_confidence_sample_1.json").write_text(
                json.dumps({"ranking_score": 0.90, "ptm": 0.88, "iptm": 0.86, "plddt": 92.0, "gpde": 0.9, "has_clash": False}),
                encoding="utf-8",
            )

            flat_structs = tmp_dir / "structs"
            out_tsv = tmp_dir / "protenix.tsv"
            out_json = tmp_dir / "protenix.json"
            subprocess.run(
                [
                    sys.executable,
                    str(PROTENIX_DIR / "bin" / "compile_protenix_outputs.py"),
                    "--out-dir",
                    str(tmp_dir / "protenix_raw"),
                    "--manifest",
                    str(manifest),
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
            self.assertEqual(rows[0]["design_id"], "bb1_mpnn_1")
            self.assertEqual(rows[0]["parent_design_id"], "bb1")
            self.assertAlmostEqual(rows[0]["confidence_score"], 0.90, places=3)
            self.assertAlmostEqual(rows[0]["complex_plddt"], 0.92, places=3)
            self.assertIn("data_s1", (flat_structs / "bb1_mpnn_1.cif").read_text(encoding="utf-8"))


class TestChai1Scripts(unittest.TestCase):
    def test_prepare_and_compile_chai1_npz(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            seq_tsv = tmp_dir / "sequences.tsv"
            seq_tsv.write_text(
                "design_id\tparent_design_id\tsequence\tparent_structure_path\n"
                "bb2_mpnn_1\tbb2\tMKTAYIAKQRQ\t/tmp/bb2.pdb\n",
                encoding="utf-8",
            )
            chai_in = tmp_dir / "chai_in"
            manifest = tmp_dir / "manifest.json"
            subprocess.run(
                [
                    sys.executable,
                    str(CHAI1_DIR / "bin" / "prepare_chai1_inputs.py"),
                    "--sequences-tsv",
                    str(seq_tsv),
                    "--target-sequence",
                    "ACDEFGHIKLMNPQRSTVWY",
                    "--ligand-smiles",
                    "CCO",
                    "--out-dir",
                    str(chai_in),
                    "--out-manifest",
                    str(manifest),
                ],
                check=True,
            )
            fa_text = (chai_in / "bb2_mpnn_1.fasta").read_text(encoding="utf-8")
            self.assertIn(">protein|name=design_chain_1", fa_text)
            self.assertIn(">protein|name=target", fa_text)
            self.assertIn(">ligand|name=ligand\nCCO", fa_text)

            job_out = tmp_dir / "chai_raw" / "bb2_mpnn_1"
            job_out.mkdir(parents=True)
            (job_out / "pred.model_idx_0.cif").write_text("data_m0\n#\n", encoding="utf-8")
            (job_out / "pred.model_idx_1.cif").write_text("data_m1\n#\n", encoding="utf-8")
            _write_chai_npz(job_out / "scores.model_idx_0.npz", agg=0.62, ptm=0.70, iptm=0.60, clash=False)
            _write_chai_npz(job_out / "scores.model_idx_1.npz", agg=0.86, ptm=0.85, iptm=0.84, clash=False)

            flat_structs = tmp_dir / "chai_structs"
            out_tsv = tmp_dir / "chai.tsv"
            out_json = tmp_dir / "chai.json"
            subprocess.run(
                [
                    sys.executable,
                    str(CHAI1_DIR / "bin" / "compile_chai1_outputs.py"),
                    "--out-dir",
                    str(tmp_dir / "chai_raw"),
                    "--manifest",
                    str(manifest),
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
            self.assertEqual(rows[0]["design_id"], "bb2_mpnn_1")
            self.assertEqual(rows[0]["parent_design_id"], "bb2")
            self.assertAlmostEqual(rows[0]["confidence_score"], 0.86, places=2)
            self.assertAlmostEqual(rows[0]["iptm"], 0.84, places=2)
            self.assertEqual(rows[0]["has_inter_chain_clashes"], False)
            self.assertIn("data_m1", (flat_structs / "bb2_mpnn_1.cif").read_text(encoding="utf-8"))


class TestProtenixChai1ContractsAndIntegration(unittest.TestCase):
    def _assert_yaml_and_params_parity(self, wrapper_dir: Path, expected_id: str) -> None:
        yaml_text = (wrapper_dir / "wrapper.yaml").read_text(encoding="utf-8")
        params_json = json.loads((wrapper_dir / "params.json").read_text(encoding="utf-8"))
        self.assertIn(f"id: {expected_id}", yaml_text)
        self.assertIn("primary: true", yaml_text)
        self.assertNotIn("\n    default:", yaml_text)

        in_params = False
        yaml_keys: set[str] = set()
        for line in yaml_text.splitlines():
            if line.startswith("params:"):
                in_params = True
                continue
            if in_params and line and not line.startswith(" "):
                break
            if in_params:
                m = re.match(r"^  ([a-z][a-z0-9_]*):\s*$", line)
                if m:
                    yaml_keys.add(m.group(1))
        self.assertEqual(yaml_keys, set(params_json.keys()), f"Param mismatch in {wrapper_dir}")

    def test_protenix_chai1_and_four_predictor_integration(self) -> None:
        self._assert_yaml_and_params_parity(PROTENIX_DIR / "wrapper", "local/modules/protenix")
        self._assert_yaml_and_params_parity(CHAI1_DIR / "wrapper", "local/modules/chai1")
        self._assert_yaml_and_params_parity(
            WRAPPERS_BASE / "subworkflows" / "local" / "seq_fold_eval" / "wrapper",
            "local/subworkflows/seq-fold-eval",
        )
        self._assert_yaml_and_params_parity(
            WRAPPERS_BASE / "workflows" / "local" / "protein_binder_design" / "wrapper",
            "local/workflows/protein-binder-design",
        )

        sfe_text = (WRAPPERS_BASE / "subworkflows" / "local" / "seq_fold_eval" / "main.nf").read_text(encoding="utf-8")
        for proc_name in ("BOLTZ2", "ALPHAFOLD3", "PROTENIX", "CHAI1"):
            self.assertIn(proc_name, sfe_text)


if __name__ == "__main__":
    unittest.main()
