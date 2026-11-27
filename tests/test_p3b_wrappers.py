import csv
import gzip
import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]


class TestP3BWrappers(unittest.TestCase):
    def test_01_vina_parser_outputs_chem_docking_eval_table(self):
        script = REPO_ROOT / "resources/wrappers/modules/local/autodock_vina/bin/parse_vina_results.py"
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            docked_dir = tmp_path / "docked"
            docked_dir.mkdir()
            (docked_dir / "ligA_out.pdbqt").write_text(
                "MODEL 1\n"
                "REMARK VINA RESULT:    -9.4      0.000      0.000\n"
                "ENDMDL\n"
                "MODEL 2\n"
                "REMARK VINA RESULT:    -8.1      1.850      3.120\n"
                "ENDMDL\n",
                encoding="utf-8",
            )
            (docked_dir / "ligB_out.pdbqt").write_text(
                "MODEL 1\n"
                "REMARK VINA RESULT:   -10.2      0.000      0.000\n"
                "ENDMDL\n",
                encoding="utf-8",
            )
            out_tsv = tmp_path / "vina.score_table.tsv"
            out_json = tmp_path / "vina.summary.json"
            subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(script),
                    "--docked_dir",
                    str(docked_dir),
                    "--output_tsv",
                    str(out_tsv),
                    "--output_json",
                    str(out_json),
                ],
                check=True,
            )
            rows = list(csv.DictReader(out_tsv.open("r", encoding="utf-8"), delimiter="\t"))
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[0]["ligand_id"], "ligB")
            self.assertAlmostEqual(float(rows[0]["score"]), -10.2)
            summary = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertEqual(summary["n_ligands"], 2)
            self.assertAlmostEqual(summary["best_affinity_kcal_mol"], -10.2)

    def test_02_gnina_sdf_parser_handles_gzipped_sdf(self):
        script = REPO_ROOT / "resources/wrappers/modules/local/gnina/bin/parse_gnina_sdf.py"
        sdf_body = (
            "cmpd_1\n"
            "  Mrv2014 01012600003D\n\n"
            "  0  0  0     0  0            999 V2000\n"
            "M  END\n"
            ">  <minimizedAffinity>\n-7.85\n\n"
            ">  <CNNscore>\n0.88\n\n"
            ">  <CNNaffinity>\n6.95\n\n"
            ">  <CNN_VS>\n6.116\n\n"
            "$$$$\n"
            "cmpd_2\n"
            "  Mrv2014 01012600003D\n\n"
            "  0  0  0     0  0            999 V2000\n"
            "M  END\n"
            ">  <minimizedAffinity>\n-8.40\n\n"
            ">  <CNNscore>\n0.92\n\n"
            ">  <CNNaffinity>\n7.40\n\n"
            "$$$$\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            sdf_gz = tmp_path / "poses.docked.sdf.gz"
            with gzip.open(sdf_gz, "wt", encoding="utf-8") as fh:
                fh.write(sdf_body)
            out_tsv = tmp_path / "gnina.score_table.tsv"
            out_json = tmp_path / "gnina.summary.json"
            subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(script),
                    "--sdf",
                    str(sdf_gz),
                    "--output_tsv",
                    str(out_tsv),
                    "--output_json",
                    str(out_json),
                ],
                check=True,
            )
            rows = list(csv.DictReader(out_tsv.open("r", encoding="utf-8"), delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["ligand_id"], "cmpd_2")
            self.assertAlmostEqual(float(rows[0]["cnn_affinity"]), 7.40)
            self.assertAlmostEqual(float(rows[0]["cnn_vs"]), 6.808)
            summary = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertEqual(summary["n_ligands"], 2)

    def test_03_fpocket_parser_computes_druggable_pocket_box(self):
        script = REPO_ROOT / "resources/wrappers/modules/local/fpocket/bin/parse_fpocket_info.py"
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            out_dir = tmp_path / "target_out"
            pockets_dir = out_dir / "pockets"
            pockets_dir.mkdir(parents=True)
            (out_dir / "target_info.txt").write_text(
                "Pocket 1 :\n"
                "\tScore : \t0.412\n"
                "\tDruggability Score : \t0.250\n"
                "\tNumber of Alpha Spheres : \t42\n"
                "\tTotal SASA : \t210.5\n"
                "\tVolume : \t450.2\n"
                "\tHydrophobicity score: \t12.4\n"
                "\tPolarity score: \t5.1\n\n"
                "Pocket 2 :\n"
                "\tScore : \t0.815\n"
                "\tDruggability Score : \t0.910\n"
                "\tNumber of Alpha Spheres : \t68\n"
                "\tTotal SASA : \t315.0\n"
                "\tVolume : \t680.4\n"
                "\tHydrophobicity score: \t28.9\n"
                "\tPolarity score: \t7.2\n",
                encoding="utf-8",
            )
            (pockets_dir / "pocket2_atm.pdb").write_text(
                "ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00  0.00           C\n"
                "ATOM      2  CA  LEU A   2      14.000  24.000  36.000  1.00  0.00           C\n",
                encoding="utf-8",
            )
            out_tsv = tmp_path / "fpocket.pockets.tsv"
            out_json = tmp_path / "fpocket.top_pocket_box.json"
            subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(script),
                    "--fpocket_out",
                    str(out_dir),
                    "--padding",
                    "6.0",
                    "--output_tsv",
                    str(out_tsv),
                    "--output_json",
                    str(out_json),
                ],
                check=True,
            )
            rows = list(csv.DictReader(out_tsv.open("r", encoding="utf-8"), delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["pocket_id"], "2")
            self.assertAlmostEqual(float(rows[0]["center_x"]), 12.0)
            self.assertAlmostEqual(float(rows[0]["center_y"]), 22.0)
            self.assertAlmostEqual(float(rows[0]["center_z"]), 33.0)
            box = json.loads(out_json.read_text(encoding="utf-8"))["top_pocket"]
            self.assertGreaterEqual(float(box["size_z"]), 18.0)

    def test_04_wrapper_yaml_and_nextflow_includes_resolve_on_disk(self):
        rel_dirs = [
            "resources/wrappers/modules/local/autodock_vina",
            "resources/wrappers/modules/local/gnina",
            "resources/wrappers/modules/local/fpocket",
            "resources/wrappers/workflows/local/rnasplice",
            "resources/wrappers/workflows/local/cutandrun",
        ]
        include_re = re.compile(r"include\s*\{([^}]+)\}\s*from\s*['\"]([^'\"]+)['\"]")
        for rel_dir in rel_dirs:
            wdir = REPO_ROOT / rel_dir
            wyaml = yaml.safe_load((wdir / "wrapper/wrapper.yaml").read_text(encoding="utf-8"))
            self.assertEqual(wyaml["schemaVersion"], 1)
            self.assertTrue((wdir / wyaml["entrypoint"]).is_file())
            self.assertTrue((wdir / wyaml["config"]).is_file())
            defaults = json.loads((wdir / wyaml["defaultsFile"]).read_text(encoding="utf-8"))
            self.assertIsInstance(defaults, dict)
            for nf in [wdir / "main.nf", wdir / wyaml["entrypoint"]]:
                txt = nf.read_text(encoding="utf-8")
                for syms_raw, inc in include_re.findall(txt):
                    resolved = (nf.parent / inc).resolve()
                    self.assertTrue(resolved.is_file(), f"Missing include {inc} from {nf}")
                    inc_txt = resolved.read_text(encoding="utf-8")
                    for sym in [s.strip().split()[-1] for s in syms_raw.split(",") if s.strip()]:
                        self.assertRegex(
                            inc_txt,
                            rf"\b(?:process|workflow)\s+{re.escape(sym)}\s*\{{",
                            f"Symbol {sym} not declared in {resolved}",
                        )


if __name__ == "__main__":
    unittest.main()
