"""Unit and contract tests for Phase 1 & Phase 2 protein design Nextflow modules:
  - resources/wrappers/modules/local/protein_metrics
  - resources/wrappers/modules/local/design_filter
  - resources/wrappers/modules/local/rfd3
  - resources/wrappers/modules/local/ligandmpnn
  - resources/wrappers/modules/local/boltz2
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
WRAPPERS_ROOT = REPO_ROOT / "resources" / "wrappers" / "modules" / "local"

PROTEIN_METRICS_DIR = WRAPPERS_ROOT / "protein_metrics"
DESIGN_FILTER_DIR = WRAPPERS_ROOT / "design_filter"
RFD3_DIR = WRAPPERS_ROOT / "rfd3"
LIGANDMPNN_DIR = WRAPPERS_ROOT / "ligandmpnn"
BOLTZ2_DIR = WRAPPERS_ROOT / "boltz2"

PROTEIN_METRICS_SCRIPT = PROTEIN_METRICS_DIR / "bin" / "compute_protein_metrics.py"
DESIGN_FILTER_SCRIPT = DESIGN_FILTER_DIR / "bin" / "merge_and_filter_designs.py"


def _format_pdb_atom(
    record: str,
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"{record:<6}{serial:>5} {atom_name:<4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 20.00          {element:>2}\n"
    )


def _write_synthetic_pdb(
    path: Path,
    shift: tuple[float, float, float] = (0.0, 0.0, 0.0),
    rotate_90_z: bool = False,
    perturb_ca: float = 0.0,
) -> None:
    """Create a 2-chain complex (chain A binder, chain B target) + 1 ligand HETATM."""
    raw_atoms = [
        ("ATOM", "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        ("ATOM", "CA", "ALA", "A", 1, 1.45, 0.0, 0.0, "C"),
        ("ATOM", "C", "ALA", "A", 1, 2.0, 1.4, 0.0, "C"),
        ("ATOM", "O", "ALA", "A", 1, 1.3, 2.4, 0.0, "O"),
        ("ATOM", "N", "GLY", "A", 2, 3.3, 1.5, 0.0, "N"),
        ("ATOM", "CA", "GLY", "A", 2, 4.0, 2.8, 0.0, "C"),
        ("ATOM", "C", "GLY", "A", 2, 5.5, 2.6, 0.0, "C"),
        ("ATOM", "O", "GLY", "A", 2, 6.0, 1.5, 0.0, "O"),
        ("ATOM", "N", "SER", "A", 3, 6.2, 3.7, 0.0, "N"),
        ("ATOM", "CA", "SER", "A", 3, 7.6, 3.7, 0.5, "C"),
        ("ATOM", "C", "SER", "A", 3, 8.2, 5.0, 0.5, "C"),
        ("ATOM", "O", "SER", "A", 3, 7.5, 6.0, 0.5, "O"),
        ("ATOM", "N", "ASP", "B", 1, 1.5, 5.2, 0.0, "N"),
        ("ATOM", "CA", "ASP", "B", 1, 2.8, 5.6, 0.0, "C"),
        ("ATOM", "C", "ASP", "B", 1, 3.8, 4.5, 0.0, "C"),
        ("ATOM", "O", "ASP", "B", 1, 3.3, 4.3, 0.0, "O"),
        ("ATOM", "N", "LYS", "B", 2, 5.1, 4.8, 0.0, "N"),
        ("ATOM", "CA", "LYS", "B", 2, 6.1, 5.9, 0.5, "C"),
        ("HETATM", "C1", "LIG", "L", 1, 4.0, 2.8, 3.2, "C"),
        ("HETATM", "O1", "LIG", "L", 1, 5.0, 2.8, 3.2, "O"),
    ]
    sx, sy, sz = shift
    lines: list[str] = []
    for idx, (rec, aname, resname, chain, resseq, x, y, z, elem) in enumerate(raw_atoms, start=1):
        if perturb_ca and chain == "A" and aname == "CA" and resseq == 3:
            z += perturb_ca
        if rotate_90_z:
            x, y = -y, x
        x += sx
        y += sy
        z += sz
        lines.append(_format_pdb_atom(rec, idx, aname, resname, chain, resseq, x, y, z, elem))
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def _write_synthetic_cif(path: Path) -> None:
    content = """data_design_cif
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.type_symbol
ATOM N ALA A 1 0.000 0.000 0.000 N
ATOM CA ALA A 1 1.450 0.000 0.000 C
ATOM C ALA A 1 2.000 1.400 0.000 C
ATOM O ALA A 1 1.300 2.400 0.000 O
ATOM N ASP B 1 1.500 5.200 0.000 N
ATOM O ASP B 1 1.300 5.400 0.000 O
"""
    path.write_text(content, encoding="utf-8")


class TestProteinMetrics(unittest.TestCase):
    def test_pdb_kabsch_rmsd_and_interface_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            ref_pdb = tmp_dir / "ref.pdb"
            rigid_pdb = tmp_dir / "design_rigid.pdb"
            perturbed_pdb = tmp_dir / "design_perturbed.pdb"
            cif_file = tmp_dir / "design_cif.cif"

            _write_synthetic_pdb(ref_pdb)
            _write_synthetic_pdb(rigid_pdb, shift=(12.5, -7.0, 3.2), rotate_90_z=True)
            _write_synthetic_pdb(perturbed_pdb, perturb_ca=3.0)
            _write_synthetic_cif(cif_file)

            out_tsv = tmp_dir / "metrics.tsv"
            out_json = tmp_dir / "metrics.json"
            proc = subprocess.run(
                [
                    sys.executable,
                    str(PROTEIN_METRICS_SCRIPT),
                    "--structures",
                    str(rigid_pdb),
                    str(perturbed_pdb),
                    str(cif_file),
                    "--reference",
                    str(ref_pdb),
                    "--binder-chains",
                    "A",
                    "--target-chains",
                    "B",
                    "--out-tsv",
                    str(out_tsv),
                    "--out-json",
                    str(out_json),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertEqual(proc.returncode, 0)

            rows = json.loads(out_json.read_text(encoding="utf-8"))
            by_id = {r["design_id"]: r for r in rows}
            self.assertEqual(set(by_id.keys()), {"design_rigid", "design_perturbed", "design_cif"})

            rigid = by_id["design_rigid"]
            self.assertAlmostEqual(rigid["ca_rmsd"], 0.0, places=3)
            self.assertAlmostEqual(rigid["backbone_rmsd"], 0.0, places=3)
            self.assertGreater(rigid["interface_contacts"], 0)
            self.assertGreater(rigid["interface_hbonds"], 0)
            self.assertGreater(rigid["ligand_contacts"], 0)
            self.assertEqual(rigid["binder_length"], 3)
            self.assertEqual(rigid["target_length"], 2)
            self.assertIsNotNone(rigid["radius_of_gyration"])

            pert = by_id["design_perturbed"]
            self.assertGreater(pert["ca_rmsd"], 0.5)

            cif_row = by_id["design_cif"]
            self.assertEqual(cif_row["num_chains"], 2)
            self.assertEqual(cif_row["binder_length"], 1)

            with out_tsv.open("r", encoding="utf-8") as fh:
                tsv_rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(len(tsv_rows), 3)


class TestDesignFilter(unittest.TestCase):
    def test_lineage_merge_dsl_three_valued_logic_and_diversity_selection(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            step1 = tmp_dir / "af3.tsv"
            step2 = tmp_dir / "metrics.tsv"

            step1.write_text(
                "design_id\tparent_design_id\tmin_pae_interaction\tplddt\n"
                "d1\tp1\t6.5\t91.2\n"
                "d2\tp1\t7.2\t88.0\n"
                "d3\tp2\t12.5\t93.0\n"
                "d4\tp2\t8.0\t86.5\n"
                "d5\tp3\t\t84.0\n",
                encoding="utf-8",
            )
            step2.write_text(
                "design_id\tca_rmsd\tinterface_contacts\tclash_count\tfallback_score\n"
                "d1\t1.1\t24\t0\t0.5\n"
                "d2\t1.4\t20\t0\t0.4\n"
                "d3\t1.0\t30\t0\t0.9\n"
                "d4\t1.8\t18\t0\t0.8\n"
                "d5\t1.5\t22\t0\t0.95\n",
                encoding="utf-8",
            )

            out_merged = tmp_dir / "merged.tsv"
            out_filtered = tmp_dir / "filtered.tsv"
            out_summary = tmp_dir / "summary.json"

            expr = (
                "((af3.min_pae_interaction <= 10.0) or (fallback_score >= 0.9)) "
                "and ca_rmsd <= 2.0 and not (clash_count > 0)"
            )
            subprocess.run(
                [
                    sys.executable,
                    str(DESIGN_FILTER_SCRIPT),
                    "--tables",
                    f"af3:{step1}",
                    str(step2),
                    "--expr",
                    expr,
                    "--rank-by",
                    "af3.plddt",
                    "--rank-direction",
                    "max",
                    "--group-by",
                    "af3.parent_design_id",
                    "--max-per-group",
                    "1",
                    "--top-n",
                    "3",
                    "--out-merged",
                    str(out_merged),
                    "--out-filtered",
                    str(out_filtered),
                    "--out-summary",
                    str(out_summary),
                ],
                capture_output=True,
                text=True,
                check=True,
            )

            summary = json.loads(out_summary.read_text(encoding="utf-8"))
            self.assertEqual(summary["total_designs"], 5)
            self.assertEqual(summary["passed_filter"], 5)
            self.assertEqual(summary["selected_designs"], 3)

            with out_filtered.open("r", encoding="utf-8") as fh:
                selected_ids = [r["design_id"] for r in csv.DictReader(fh, delimiter="\t")]
            self.assertEqual(selected_ids, ["d3", "d1", "d5"])


class TestRfd3Scripts(unittest.TestCase):
    def test_prepare_and_compile_rfd3(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            target_pdb = tmp_dir / "target.pdb"
            _write_synthetic_pdb(target_pdb)

            staged_json = tmp_dir / "staged.json"
            subprocess.run(
                [
                    sys.executable,
                    str(RFD3_DIR / "bin" / "prepare_rfd3_inputs.py"),
                    "--design-name",
                    "binder",
                    "--contig",
                    "A1-100/0 60-80",
                    "--input-structure",
                    str(target_pdb),
                    "--ligand",
                    "LIG",
                    "--chunk-id",
                    "2",
                    "--out-json",
                    str(staged_json),
                ],
                check=True,
            )
            spec = json.loads(staged_json.read_text(encoding="utf-8"))
            self.assertIn("binder_chunk2", spec)
            self.assertEqual(spec["binder_chunk2"]["contig"], "A1-100/0 60-80")
            self.assertEqual(spec["binder_chunk2"]["ligand"], "LIG")

            raw_dir = tmp_dir / "rfd3_raw"
            raw_dir.mkdir()
            _write_synthetic_cif(raw_dir / "binder_chunk2_0_model_0.cif")
            (raw_dir / "binder_chunk2_0_model_0.json").write_text(
                json.dumps({"diffused_length": 72, "sampler": {"timesteps": 200}}),
                encoding="utf-8",
            )

            out_tsv = tmp_dir / "rfd3.tsv"
            out_json = tmp_dir / "rfd3.json"
            subprocess.run(
                [
                    sys.executable,
                    str(RFD3_DIR / "bin" / "compile_rfd3_outputs.py"),
                    "--out-dir",
                    str(raw_dir),
                    "--out-tsv",
                    str(out_tsv),
                    "--out-json",
                    str(out_json),
                ],
                check=True,
            )
            rows = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["design_id"], "binder_chunk2_0_model_0")
            self.assertEqual(rows[0]["diffused_length"], 72)
            self.assertEqual(rows[0]["sampler.timesteps"], 200)


class TestLigandMpnnAndBoltz2Bridge(unittest.TestCase):
    def test_ligandmpnn_to_boltz2_prepare_and_compile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            backbone = tmp_dir / "bb_0.pdb"
            _write_synthetic_pdb(backbone)

            multi_json = tmp_dir / "multi.json"
            mpnn_manifest = tmp_dir / "mpnn_manifest.json"
            subprocess.run(
                [
                    sys.executable,
                    str(LIGANDMPNN_DIR / "bin" / "prepare_ligandmpnn_inputs.py"),
                    "--structures",
                    str(backbone),
                    "--out-multi-json",
                    str(multi_json),
                    "--out-manifest",
                    str(mpnn_manifest),
                ],
                check=True,
            )
            self.assertIn(str(backbone.resolve()), json.loads(multi_json.read_text(encoding="utf-8")))

            mpnn_raw = tmp_dir / "mpnn_raw"
            (mpnn_raw / "seqs").mkdir(parents=True)
            (mpnn_raw / "seqs" / "bb_0.fa").write_text(
                ">bb_0, T=0.1, seed=111\n"
                "MKTAYIAKQRQ\n"
                ">bb_0, id=1, T=0.1, seed=111, overall_confidence=0.48, ligand_confidence=0.52, seq_rec=0.45\n"
                "MKTAYIAKQRQ\n"
                ">bb_0, id=2, T=0.1, seed=111, overall_confidence=0.41, ligand_confidence=0.44, seq_rec=0.39\n"
                "MKTGYVAKQRQ\n",
                encoding="utf-8",
            )

            seq_tsv = tmp_dir / "sequences.tsv"
            seq_fa = tmp_dir / "sequences.fasta"
            seq_json = tmp_dir / "mpnn_metrics.json"
            subprocess.run(
                [
                    sys.executable,
                    str(LIGANDMPNN_DIR / "bin" / "compile_ligandmpnn_outputs.py"),
                    "--out-dir",
                    str(mpnn_raw),
                    "--manifest",
                    str(mpnn_manifest),
                    "--out-tsv",
                    str(seq_tsv),
                    "--out-fasta",
                    str(seq_fa),
                    "--out-json",
                    str(seq_json),
                ],
                check=True,
            )
            mpnn_rows = json.loads(seq_json.read_text(encoding="utf-8"))
            self.assertEqual(len(mpnn_rows), 2)
            self.assertEqual(mpnn_rows[0]["design_id"], "bb_0_mpnn_1")
            self.assertEqual(mpnn_rows[0]["parent_design_id"], "bb_0")
            self.assertAlmostEqual(mpnn_rows[0]["ligand_confidence"], 0.52, places=3)

            boltz_in = tmp_dir / "boltz_inputs"
            boltz_manifest = tmp_dir / "boltz_manifest.json"
            subprocess.run(
                [
                    sys.executable,
                    str(BOLTZ2_DIR / "bin" / "prepare_boltz2_inputs.py"),
                    "--sequences-tsv",
                    str(seq_tsv),
                    "--target-sequence",
                    "ACDEFGHIKLMNPQRSTVWY",
                    "--ligand-ccd",
                    "ATP",
                    "--predict-affinity",
                    "--out-dir",
                    str(boltz_in),
                    "--out-manifest",
                    str(boltz_manifest),
                ],
                check=True,
            )
            yaml1 = (boltz_in / "bb_0_mpnn_1.yaml").read_text(encoding="utf-8")
            self.assertIn("ccd: ATP", yaml1)
            self.assertIn("affinity:", yaml1)
            self.assertIn("ACDEFGHIKLMNPQRSTVWY", yaml1)

            job_dir = tmp_dir / "boltz_raw" / "boltz_results_inputs" / "predictions" / "bb_0_mpnn_1"
            job_dir.mkdir(parents=True)
            _write_synthetic_cif(job_dir / "bb_0_mpnn_1_model_0.cif")
            _write_synthetic_cif(job_dir / "bb_0_mpnn_1_model_1.cif")
            (job_dir / "confidence_bb_0_mpnn_1_model_0.json").write_text(
                json.dumps({"confidence_score": 0.72, "iptm": 0.65, "complex_plddt": 0.80}),
                encoding="utf-8",
            )
            (job_dir / "confidence_bb_0_mpnn_1_model_1.json").write_text(
                json.dumps({"confidence_score": 0.89, "iptm": 0.84, "complex_plddt": 0.91}),
                encoding="utf-8",
            )
            (job_dir / "affinity_bb_0_mpnn_1.json").write_text(
                json.dumps({"affinity_pred_value": -2.1, "affinity_probability_binary": 0.93}),
                encoding="utf-8",
            )

            flat_structs = tmp_dir / "boltz_structs"
            boltz_tsv = tmp_dir / "boltz.tsv"
            boltz_json = tmp_dir / "boltz.json"
            subprocess.run(
                [
                    sys.executable,
                    str(BOLTZ2_DIR / "bin" / "compile_boltz2_outputs.py"),
                    "--out-dir",
                    str(tmp_dir / "boltz_raw"),
                    "--manifest",
                    str(boltz_manifest),
                    "--out-structures-dir",
                    str(flat_structs),
                    "--out-tsv",
                    str(boltz_tsv),
                    "--out-json",
                    str(boltz_json),
                ],
                check=True,
            )
            b_rows = json.loads(boltz_json.read_text(encoding="utf-8"))
            self.assertEqual(len(b_rows), 1)
            self.assertEqual(b_rows[0]["design_id"], "bb_0_mpnn_1")
            self.assertEqual(b_rows[0]["parent_design_id"], "bb_0")
            self.assertAlmostEqual(b_rows[0]["confidence_score"], 0.89, places=3)
            self.assertAlmostEqual(b_rows[0]["affinity_probability_binary"], 0.93, places=3)
            self.assertTrue((flat_structs / "bb_0_mpnn_1.cif").is_file())


class TestWrapperContractCompliance(unittest.TestCase):
    def test_all_protein_design_modules_satisfy_wrapper_contract_files(self) -> None:
        expected = [
            (PROTEIN_METRICS_DIR, "local/modules/protein-metrics", ["compute_protein_metrics.py"]),
            (DESIGN_FILTER_DIR, "local/modules/design-filter", ["merge_and_filter_designs.py"]),
            (RFD3_DIR, "local/modules/rfd3", ["prepare_rfd3_inputs.py", "compile_rfd3_outputs.py"]),
            (LIGANDMPNN_DIR, "local/modules/ligandmpnn", ["prepare_ligandmpnn_inputs.py", "compile_ligandmpnn_outputs.py"]),
            (BOLTZ2_DIR, "local/modules/boltz2", ["prepare_boltz2_inputs.py", "compile_boltz2_outputs.py"]),
        ]
        for mod_dir, expected_id, scripts in expected:
            self.assertTrue((mod_dir / "main.nf").is_file(), f"missing main.nf in {mod_dir}")
            self.assertTrue((mod_dir / "meta.yml").is_file(), f"missing meta.yml in {mod_dir}")
            self.assertTrue((mod_dir / "environment.yml").is_file(), f"missing environment.yml in {mod_dir}")
            for script_name in scripts:
                self.assertTrue((mod_dir / "bin" / script_name).is_file(), f"missing bin/{script_name}")

            wrapper_dir = mod_dir / "wrapper"
            yaml_text = (wrapper_dir / "wrapper.yaml").read_text(encoding="utf-8")
            params_json = json.loads((wrapper_dir / "params.json").read_text(encoding="utf-8"))
            self.assertTrue((wrapper_dir / "main.nf").is_file())
            self.assertTrue((wrapper_dir / "nextflow.config").is_file())

            self.assertIn(f"id: {expected_id}", yaml_text)
            self.assertIn("primary: true", yaml_text)
            self.assertNotIn("\n    default:", yaml_text, "defaults must live only in wrapper/params.json")
            self.assertIn("outdir", params_json)
            self.assertIsInstance(params_json, dict)


if __name__ == "__main__":
    unittest.main()
