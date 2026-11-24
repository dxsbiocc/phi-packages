"""Unit and contract tests for extended protein design Nextflow modules:
  - resources/wrappers/modules/local/bindcraft
  - resources/wrappers/modules/local/foldseek
  - resources/wrappers/modules/local/esmfold2
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
LOCAL_MODS = REPO_ROOT / "resources" / "wrappers" / "modules" / "local"

BINDCRAFT_DIR = LOCAL_MODS / "bindcraft"
FOLDSEEK_DIR = LOCAL_MODS / "foldseek"
ESMFOLD2_DIR = LOCAL_MODS / "esmfold2"


class TestBindCraftScripts(unittest.TestCase):
    def test_prepare_and_compile_bindcraft(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            target_pdb = tmp_dir / "target.pdb"
            target_pdb.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 85.00           C\nEND\n", encoding="utf-8")

            design_dir = tmp_dir / "bindcraft_raw"
            staged_json = tmp_dir / "settings.json"
            subprocess.run(
                [
                    sys.executable,
                    str(BINDCRAFT_DIR / "bin" / "prepare_bindcraft_inputs.py"),
                    "--starting-pdb",
                    str(target_pdb),
                    "--design-path",
                    str(design_dir),
                    "--binder-name",
                    "pdl1_binder",
                    "--chains",
                    "A",
                    "--target-hotspot-residues",
                    "A56,A58-62",
                    "--min-length",
                    "70",
                    "--max-length",
                    "110",
                    "--num-designs",
                    "5",
                    "--out-json",
                    str(staged_json),
                ],
                check=True,
            )
            cfg = json.loads(staged_json.read_text(encoding="utf-8"))
            self.assertEqual(cfg["binder_name"], "pdl1_binder")
            self.assertEqual(cfg["lengths"], [70, 110])
            self.assertTrue(Path(cfg["starting_pdb"]).is_absolute())

            # Simulate BindCraft outputs
            (design_dir / "Accepted").mkdir(parents=True)
            (design_dir / "Accepted" / "pdl1_binder_l85_s123_mpnn1.pdb").write_text("END\n", encoding="utf-8")
            (design_dir / "final_design_stats.csv").write_text(
                "Design,Sequence,Average_pLDDT,Average_i_pTM,Average_pTM,Average_i_pAE,Average_dG,Average_ShapeComplementarity,Average_n_InterfaceResidues,Average_n_InterfaceHbonds,Average_Relaxed_Clashes\n"
                "pdl1_binder_l85_s123_mpnn1,MKTAYIAKQRQ,91.5,0.84,0.88,0.18,-42.3,0.72,22,5,0\n",
                encoding="utf-8",
            )

            flat_structs = tmp_dir / "bc_structs"
            out_tsv = tmp_dir / "bc.tsv"
            out_json = tmp_dir / "bc.json"
            subprocess.run(
                [
                    sys.executable,
                    str(BINDCRAFT_DIR / "bin" / "compile_bindcraft_outputs.py"),
                    "--design-path",
                    str(design_dir),
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
            self.assertEqual(rows[0]["design_id"], "pdl1_binder_l85_s123_mpnn1")
            self.assertEqual(rows[0]["parent_design_id"], "pdl1_binder_l85_s123")
            self.assertAlmostEqual(rows[0]["complex_plddt"], 0.915, places=3)
            self.assertAlmostEqual(rows[0]["iptm"], 0.84, places=3)
            self.assertTrue((flat_structs / "pdl1_binder_l85_s123_mpnn1.pdb").is_file())


class TestFoldseekScripts(unittest.TestCase):
    def test_compile_search_and_cluster(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            qdir = tmp_dir / "queries"
            qdir.mkdir()
            (qdir / "d1.pdb").write_text("END\n", encoding="utf-8")
            (qdir / "d2.pdb").write_text("END\n", encoding="utf-8")
            (qdir / "d3.pdb").write_text("END\n", encoding="utf-8")

            # Search mode: d1 has high TM hit (not novel), d2 has low TM hit (novel), d3 has 0 hits (novel)
            raw_search = tmp_dir / "hits.tsv"
            raw_search.write_text(
                "d1.pdb\t1abc_A\t0.45\t80\t1e-12\t180\t0.78\t0.75\t0.70\t0.81\n"
                "d2.pdb\t2xyz_A\t0.20\t60\t1e-2\t42\t0.38\t0.35\t0.30\t0.44\n",
                encoding="utf-8",
            )
            s_tsv = tmp_dir / "search.tsv"
            s_json = tmp_dir / "search.json"
            subprocess.run(
                [
                    sys.executable,
                    str(FOLDSEEK_DIR / "bin" / "compile_foldseek_outputs.py"),
                    "--mode",
                    "search",
                    "--input-tsv",
                    str(raw_search),
                    "--query-dir",
                    str(qdir),
                    "--novelty-tm-cutoff",
                    "0.5",
                    "--out-tsv",
                    str(s_tsv),
                    "--out-json",
                    str(s_json),
                ],
                check=True,
            )
            s_rows = {r["design_id"]: r for r in json.loads(s_json.read_text(encoding="utf-8"))}
            self.assertEqual(set(s_rows.keys()), {"d1", "d2", "d3"})
            self.assertEqual(s_rows["d1"]["is_novel"], "false")
            self.assertEqual(s_rows["d2"]["is_novel"], "true")
            self.assertEqual(s_rows["d3"]["num_hits"], 0)
            self.assertEqual(s_rows["d3"]["is_novel"], "true")

            # Cluster mode: d1+d2 in cluster 1, d3 singleton in cluster 2
            raw_clu = tmp_dir / "clu.tsv"
            raw_clu.write_text("d1.pdb\td1.pdb\nd1.pdb\td2.pdb\nd3.pdb\td3.pdb\n", encoding="utf-8")
            c_tsv = tmp_dir / "cluster.tsv"
            c_json = tmp_dir / "cluster.json"
            subprocess.run(
                [
                    sys.executable,
                    str(FOLDSEEK_DIR / "bin" / "compile_foldseek_outputs.py"),
                    "--mode",
                    "cluster",
                    "--input-tsv",
                    str(raw_clu),
                    "--out-tsv",
                    str(c_tsv),
                    "--out-json",
                    str(c_json),
                ],
                check=True,
            )
            c_rows = {r["design_id"]: r for r in json.loads(c_json.read_text(encoding="utf-8"))}
            self.assertEqual(c_rows["d1"]["cluster_id"], "clu_1")
            self.assertEqual(c_rows["d2"]["cluster_id"], "clu_1")
            self.assertEqual(c_rows["d1"]["is_representative"], "true")
            self.assertEqual(c_rows["d2"]["is_representative"], "false")
            self.assertEqual(c_rows["d1"]["cluster_size"], 2)
            self.assertEqual(c_rows["d3"]["cluster_id"], "clu_2")


class TestESMFold2Scripts(unittest.TestCase):
    def test_prepare_and_compile_esmfold2(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            seq_tsv = tmp_dir / "sequences.tsv"
            seq_tsv.write_text(
                "design_id\tparent_design_id\tsequence\tparent_structure_path\n"
                "bb0_mpnn_1\tbb0\tMKTAYIAKQRQ\t/tmp/bb0.pdb\n",
                encoding="utf-8",
            )
            staged_fa = tmp_dir / "staged.fasta"
            manifest = tmp_dir / "manifest.json"
            subprocess.run(
                [
                    sys.executable,
                    str(ESMFOLD2_DIR / "bin" / "prepare_esmfold2_inputs.py"),
                    "--sequences-tsv",
                    str(seq_tsv),
                    "--out-fasta",
                    str(staged_fa),
                    "--out-manifest",
                    str(manifest),
                ],
                check=True,
            )
            self.assertIn(">bb0_mpnn_1", staged_fa.read_text(encoding="utf-8"))

            # Simulate ESMFold2 output PDB with pLDDT in B-factor column
            raw_dir = tmp_dir / "esm_raw"
            raw_dir.mkdir()
            (raw_dir / "bb0_mpnn_1.pdb").write_text(
                "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 84.00           C\n"
                "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00 88.00           C\n"
                "END\n",
                encoding="utf-8",
            )
            flat_structs = tmp_dir / "esm_structs"
            out_tsv = tmp_dir / "esm.tsv"
            out_json = tmp_dir / "esm.json"
            subprocess.run(
                [
                    sys.executable,
                    str(ESMFOLD2_DIR / "bin" / "compile_esmfold2_outputs.py"),
                    "--out-dir",
                    str(raw_dir),
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
            self.assertEqual(rows[0]["design_id"], "bb0_mpnn_1")
            self.assertEqual(rows[0]["parent_design_id"], "bb0")
            self.assertAlmostEqual(rows[0]["plddt"], 86.0, places=3)
            self.assertAlmostEqual(rows[0]["complex_plddt"], 0.86, places=3)
            self.assertTrue((flat_structs / "bb0_mpnn_1.pdb").is_file())


class TestExtraModulesContractCompliance(unittest.TestCase):
    def test_extra_modules_contract_and_param_parity(self) -> None:
        expected = [
            (BINDCRAFT_DIR, "local/modules/bindcraft", ["prepare_bindcraft_inputs.py", "compile_bindcraft_outputs.py"]),
            (FOLDSEEK_DIR, "local/modules/foldseek", ["compile_foldseek_outputs.py"]),
            (ESMFOLD2_DIR, "local/modules/esmfold2", ["prepare_esmfold2_inputs.py", "compile_esmfold2_outputs.py"]),
        ]
        for mod_dir, expected_id, scripts in expected:
            self.assertTrue((mod_dir / "main.nf").is_file())
            self.assertTrue((mod_dir / "meta.yml").is_file())
            self.assertTrue((mod_dir / "environment.yml").is_file())
            for s in scripts:
                self.assertTrue((mod_dir / "bin" / s).is_file())

            wdir = mod_dir / "wrapper"
            yaml_text = (wdir / "wrapper.yaml").read_text(encoding="utf-8")
            params_json = json.loads((wdir / "params.json").read_text(encoding="utf-8"))
            self.assertTrue((wdir / "main.nf").is_file())
            self.assertTrue((wdir / "nextflow.config").is_file())
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
            self.assertEqual(yaml_keys, set(params_json.keys()), f"Param mismatch in {wdir}")


if __name__ == "__main__":
    unittest.main()
