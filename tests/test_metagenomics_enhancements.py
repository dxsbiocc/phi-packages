"""Unit, contract, syntax, dependency-boundary, and end-to-end tests for
metagenomics plugin 0.2.0 enhancements:
  - mag-workflow: mg_mag_traits (AMRFinderPlus, ABRicate AMR/VFDB, RGI, antiSMASH)
  - profiling-workflow: mg_diversity (PERMDISP + multi-covariate PERMANOVA)
  - profiling-workflow: mg_network (Spearman co-occurrence network + hubs)
  - profiling-workflow: mg_classify (leak-free CV mRMR + L2 logistic regression)
  - amplicon-workflow: mg_amplicon_unifrac (Weighted & Unweighted UniFrac + PCoA/PERMANOVA/PERMDISP)
"""

from __future__ import annotations

import ast
import importlib.util
import json
import py_compile
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
MG_PLUGIN = REPO_ROOT / "resources" / "plugins" / "metagenomics"
MAG_SCRIPTS = MG_PLUGIN / "skills" / "mag-workflow" / "scripts"
PROF_SCRIPTS = MG_PLUGIN / "skills" / "profiling-workflow" / "scripts"
AMP_SCRIPTS = MG_PLUGIN / "skills" / "amplicon-workflow" / "scripts"

ALLOWED_THIRD_PARTY = {"numpy", "pandas", "scipy", "matplotlib", "_common"}
HAS_SCIPY_STACK = all(
    importlib.util.find_spec(pkg) is not None
    for pkg in ("numpy", "pandas", "scipy", "matplotlib")
)


def _run_json(cmd: list[str]) -> dict:
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return json.loads(proc.stdout)


class TestMetagenomicsStaticAndContracts(unittest.TestCase):
    def test_scripts_compile_and_respect_locked_environment(self) -> None:
        scripts = [
            MAG_SCRIPTS / "mag_traits.py",
            MAG_SCRIPTS / "mag_metabolism.py",
            PROF_SCRIPTS / "diversity.py",
            PROF_SCRIPTS / "network.py",
            PROF_SCRIPTS / "classify.py",
            PROF_SCRIPTS / "multivariable.py",
            AMP_SCRIPTS / "amplicon_unifrac.py",
        ]
        stdlib_names = set(sys.stdlib_module_names)
        for script in scripts:
            self.assertTrue(script.is_file(), f"missing script: {script}")
            py_compile.compile(str(script), doraise=True)
            tree = ast.parse(script.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        top = alias.name.split(".")[0]
                        self.assertTrue(
                            top in stdlib_names or top in ALLOWED_THIRD_PARTY,
                            f"disallowed import {top!r} in {script.name}",
                        )
                elif isinstance(node, ast.ImportFrom) and node.module:
                    top = node.module.split(".")[0]
                    self.assertTrue(
                        top in stdlib_names or top in ALLOWED_THIRD_PARTY,
                        f"disallowed import {top!r} in {script.name}",
                    )

    def test_schemas_and_plugin_manifest_wiring(self) -> None:
        manifest = (MG_PLUGIN / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("version: 0.3.0", manifest)
        self.assertIn("skills/amplicon-workflow", manifest)

        schemas = [
            MG_PLUGIN / "skills" / "mag-workflow" / "schemas" / "mag-traits-result.json",
            MG_PLUGIN / "skills" / "mag-workflow" / "schemas" / "mag-metabolism-result.json",
            MG_PLUGIN / "skills" / "profiling-workflow" / "schemas" / "diversity-result.json",
            MG_PLUGIN / "skills" / "profiling-workflow" / "schemas" / "network-result.json",
            MG_PLUGIN / "skills" / "profiling-workflow" / "schemas" / "classify-result.json",
            MG_PLUGIN / "skills" / "profiling-workflow" / "schemas" / "multivariable-result.json",
            MG_PLUGIN / "skills" / "amplicon-workflow" / "schemas" / "amplicon-unifrac-result.json",
        ]
        for sp in schemas:
            self.assertTrue(sp.is_file(), f"missing schema: {sp}")
            data = json.loads(sp.read_text(encoding="utf-8"))
            self.assertEqual(data.get("type"), "object")
            self.assertIn("required", data)
            self.assertIn("properties", data)


@unittest.skipUnless(
    HAS_SCIPY_STACK,
    "runtime execution tests require plugin:metagenomics env (numpy/pandas/scipy/matplotlib)",
)
class TestMetagenomicsRuntimeExecution(unittest.TestCase):
    def test_mg_mag_traits_amr_vfdb_rgi_antismash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            amr_tsv = tmp_dir / "bin_1.amrfinder.tsv"
            amr_tsv.write_text(
                "Gene symbol\tElement type\tClass\t% Coverage of reference sequence\t% Identity to reference sequence\n"
                "blaTEM-1\tAMR\tBETA-LACTAM\t99.5\t98.2\n"
                "iroN\tVIRULENCE\tSIDEROPHORE\t95.0\t92.0\n"
                "arsB\tSTRESS\tARSENIC\t90.0\t88.0\n",
                encoding="utf-8",
            )
            abr_tsv = tmp_dir / "abricate_vfdb.tsv"
            abr_tsv.write_text(
                "#FILE\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tRESISTANCE\n"
                "bin_2.fa\tvirB\t96.0\t91.5\tvfdb\t\n"
                "bin_2.fa\ttetA\t99.0\t97.0\tncbi\tTETRACYCLINE\n",
                encoding="utf-8",
            )
            rgi_tsv = tmp_dir / "bin_1.rgi.txt"
            rgi_tsv.write_text(
                "ORF_ID\tBest_Hit_ARO\tDrug Class\tBest_Identities\tPercentage Length of Reference Sequence\n"
                "orf_1\tmexB\tfluoroquinolone antibiotic; macrolide antibiotic\t94.0\t98.0\n",
                encoding="utf-8",
            )
            anti_dir = tmp_dir / "antismash" / "bin_2"
            anti_dir.mkdir(parents=True)
            (anti_dir / "bin_2.json").write_text(
                json.dumps(
                    {
                        "records": [
                            {
                                "features": [
                                    {"type": "region", "qualifiers": {"product": ["NRPS", "T1PKS"]}},
                                    {"type": "region", "qualifiers": {"product": ["terpene"]}},
                                ]
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )

            out_dir = tmp_dir / "out_traits"
            res = _run_json(
                [
                    sys.executable,
                    str(MAG_SCRIPTS / "mag_traits.py"),
                    "--amrfinder",
                    str(amr_tsv),
                    "--abricate",
                    str(abr_tsv),
                    "--rgi",
                    str(rgi_tsv),
                    "--antismash",
                    str(tmp_dir / "antismash"),
                    "--output_dir",
                    str(out_dir),
                ]
            )
            self.assertEqual(res["tool"], "mg_mag_traits")
            self.assertEqual(res["n_bins"], 2)
            self.assertEqual(res["total_amr_hits"], 3)
            self.assertEqual(res["total_virulence_hits"], 2)
            self.assertEqual(res["total_stress_hits"], 1)
            self.assertEqual(res["total_bgc_regions"], 2)
            self.assertTrue((out_dir / "mag-traits-summary.tsv").is_file())
            self.assertTrue((out_dir / "mag-amr-matrix.tsv").is_file())
            self.assertTrue((out_dir / "mag-bgc-matrix.tsv").is_file())
            self.assertTrue((out_dir / "mag-traits-overview.png").is_file())

    def test_mg_diversity_permdisp_and_covariate_permanova(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            rng = random.Random(42)
            samples = [f"S{i}" for i in range(12)]
            rows = {"SpA": [], "SpB": [], "SpC": [], "SpD": []}
            for i in range(12):
                if i < 6:
                    rows["SpA"].append(80 + rng.randint(0, 9))
                    rows["SpB"].append(15 + rng.randint(0, 4))
                    rows["SpC"].append(5)
                    rows["SpD"].append(2)
                else:
                    rows["SpA"].append(5)
                    rows["SpB"].append(10 + rng.randint(0, 4))
                    rows["SpC"].append(75 + rng.randint(0, 9))
                    rows["SpD"].append(10)
            mat_path = tmp_dir / "matrix.tsv"
            lines = ["taxon\t" + "\t".join(samples)]
            for taxon, vals in rows.items():
                lines.append(taxon + "\t" + "\t".join(str(v) for v in vals))
            mat_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

            ages = [35, 42, 51, 60, 48, 55, 36, 44, 52, 59, 47, 54]
            meta_lines = ["sample\tgroup\tsex\tage"]
            for i, s in enumerate(samples):
                grp = "case" if i < 6 else "ctrl"
                sex = "M" if i % 2 == 0 else "F"
                meta_lines.append(f"{s}\t{grp}\t{sex}\t{ages[i]}")
            meta_path = tmp_dir / "meta.tsv"
            meta_path.write_text("\n".join(meta_lines) + "\n", encoding="utf-8")

            out_dir = tmp_dir / "out_div"
            res = _run_json(
                [
                    sys.executable,
                    str(PROF_SCRIPTS / "diversity.py"),
                    "--input",
                    str(mat_path),
                    "--values",
                    "counts",
                    "--metadata",
                    str(meta_path),
                    "--group_col",
                    "group",
                    "--covariates",
                    "group",
                    "sex",
                    "age",
                    "--permutations",
                    "199",
                    "--output_dir",
                    str(out_dir),
                ]
            )
            self.assertIsNotNone(res["permanova"])
            self.assertLess(res["permanova"]["p_value"], 0.05)
            self.assertIsNotNone(res["permdisp"])
            self.assertEqual(len(res["covariate_permanova"]), 3)
            cov_by_name = {c["covariate"]: c for c in res["covariate_permanova"]}
            self.assertEqual(cov_by_name["age"]["kind"], "numeric")
            self.assertLess(cov_by_name["group"]["p_value"], 0.05)
            self.assertTrue((out_dir / "covariate-permanova.tsv").is_file())

    def test_mg_network_spearman_cooccurrence_and_hubs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            rng = random.Random(7)
            samples = [f"S{i}" for i in range(16)]
            base = [1.0 + i * 1.2 for i in range(16)]
            feats = {
                "Hub_A": [max(0.1, b + rng.gauss(0, 0.2)) for b in base],
                "Co_B": [max(0.1, b * 1.5 + rng.gauss(0, 0.3)) for b in base],
                "Co_C": [max(0.1, b * 0.8 + rng.gauss(0, 0.2)) for b in base],
                "Anti_D": [max(0.1, 25.0 - b + rng.gauss(0, 0.3)) for b in base],
                "Noise_E": [rng.uniform(1.0, 5.0) for _ in base],
            }
            mat_path = tmp_dir / "mat.tsv"
            lines = ["feature\t" + "\t".join(samples)]
            for f, vals in feats.items():
                lines.append(f + "\t" + "\t".join(f"{v:.4f}" for v in vals))
            mat_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

            out_dir = tmp_dir / "out_net"
            res = _run_json(
                [
                    sys.executable,
                    str(PROF_SCRIPTS / "network.py"),
                    "--input",
                    str(mat_path),
                    "--min_abs_rho",
                    "0.6",
                    "--fdr",
                    "0.05",
                    "--output_dir",
                    str(out_dir),
                ]
            )
            self.assertEqual(res["tool"], "mg_network")
            self.assertGreaterEqual(res["n_edges"], 4)
            self.assertGreaterEqual(res["n_positive_edges"], 3)
            self.assertGreaterEqual(res["n_negative_edges"], 1)
            self.assertTrue((out_dir / "cooccurrence-network.png").is_file())

    def test_mg_classify_nested_cv_mrmr_logistic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            rng = random.Random(99)
            n_per = 15
            samples = [f"C{i}" for i in range(n_per)] + [f"H{i}" for i in range(n_per)]
            lines = ["feature\t" + "\t".join(samples)]
            for f_idx in range(15):
                if f_idx < 3:
                    vals = [max(0.0, rng.gauss(12.0, 1.5)) for _ in range(n_per)] + [
                        max(0.0, rng.gauss(2.0, 1.5)) for _ in range(n_per)
                    ]
                else:
                    vals = [max(0.0, rng.gauss(5.0, 2.0)) for _ in range(2 * n_per)]
                lines.append(f"Feat_{f_idx}\t" + "\t".join(f"{v:.4f}" for v in vals))
            mat_path = tmp_dir / "feat.tsv"
            mat_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

            meta_lines = ["sample\tgroup"] + [
                f"{s}\t{'CRC' if i < n_per else 'Control'}" for i, s in enumerate(samples)
            ]
            meta_path = tmp_dir / "meta.tsv"
            meta_path.write_text("\n".join(meta_lines) + "\n", encoding="utf-8")

            out_dir = tmp_dir / "out_cls"
            res = _run_json(
                [
                    sys.executable,
                    str(PROF_SCRIPTS / "classify.py"),
                    "--input",
                    str(mat_path),
                    "--metadata",
                    str(meta_path),
                    "--case",
                    "CRC",
                    "--control",
                    "Control",
                    "--top_features",
                    "5",
                    "--folds",
                    "5",
                    "--output_dir",
                    str(out_dir),
                ]
            )
            self.assertEqual(res["tool"], "mg_classify")
            self.assertGreater(res["oof_auroc"], 0.90)
            self.assertTrue((out_dir / "roc-curve.png").is_file())

    def test_mg_amplicon_unifrac_weighted_and_unweighted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            tree_path = tmp_dir / "rooted.nwk"
            tree_path.write_text("((ASV1:0.2,ASV2:0.2):0.5,(ASV3:0.2,ASV4:0.2):0.5);", encoding="utf-8")
            asv_path = tmp_dir / "asv.tsv"
            asv_path.write_text(
                "ASV\tA1\tA2\tA3\tA4\tB1\tB2\tB3\tB4\n"
                "ASV1\t90\t80\t85\t88\t0\t0\t0\t0\n"
                "ASV2\t10\t20\t15\t12\t0\t0\t0\t0\n"
                "ASV3\t0\t0\t0\t0\t90\t75\t80\t85\n"
                "ASV4\t0\t0\t0\t0\t10\t25\t20\t15\n",
                encoding="utf-8",
            )
            meta_path = tmp_dir / "meta.tsv"
            meta_path.write_text(
                "sample\tgroup\nA1\tR\nA2\tR\nA3\tR\nA4\tR\nB1\tS\nB2\tS\nB3\tS\nB4\tS\n",
                encoding="utf-8",
            )
            out_dir = tmp_dir / "out_uni"
            res = _run_json(
                [
                    sys.executable,
                    str(AMP_SCRIPTS / "amplicon_unifrac.py"),
                    "--asv_table",
                    str(asv_path),
                    "--tree",
                    str(tree_path),
                    "--metric",
                    "weighted",
                    "--metadata",
                    str(meta_path),
                    "--permutations",
                    "199",
                    "--output_dir",
                    str(out_dir),
                ]
            )
            self.assertEqual(res["tool"], "mg_amplicon_unifrac")
            self.assertEqual(res["n_asvs_matched"], 4)
            self.assertLess(res["permanova"]["p_value"], 0.05)
            self.assertTrue((out_dir / "weighted-unifrac.tsv").is_file())
            self.assertTrue((out_dir / "unweighted-unifrac.tsv").is_file())
            self.assertTrue((out_dir / "pcoa-weighted-unifrac.png").is_file())

    def test_mg_multivariable_clr_regression(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            rng = random.Random(11)
            samples = [f"S{i}" for i in range(16)]
            lines = ["feature\t" + "\t".join(samples)]
            # Sp_pos strongly enriched in case; Sp_bg flat; Sp_age driven by age
            ages = [30 + i * 2 for i in range(16)]
            lines.append("Sp_pos\t" + "\t".join(f"{50.0 + rng.gauss(0, 2) if i < 8 else 5.0 + rng.gauss(0, 1):.3f}" for i in range(16)))
            lines.append("Sp_bg1\t" + "\t".join(f"{20.0 + rng.gauss(0, 1):.3f}" for _ in range(16)))
            lines.append("Sp_bg2\t" + "\t".join(f"{15.0 + rng.gauss(0, 1):.3f}" for _ in range(16)))
            lines.append("Sp_age\t" + "\t".join(f"{ages[i] * 0.8 + rng.gauss(0, 1):.3f}" for i in range(16)))
            mat_path = tmp_dir / "mat.tsv"
            mat_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

            meta_lines = ["sample\tgroup\tsex\tage"] + [
                f"{s}\t{'case' if i < 8 else 'ctrl'}\t{'M' if i % 2 == 0 else 'F'}\t{ages[i]}"
                for i, s in enumerate(samples)
            ]
            meta_path = tmp_dir / "meta.tsv"
            meta_path.write_text("\n".join(meta_lines) + "\n", encoding="utf-8")

            out_dir = tmp_dir / "out_mv"
            res = _run_json(
                [
                    sys.executable,
                    str(PROF_SCRIPTS / "multivariable.py"),
                    "--input",
                    str(mat_path),
                    "--metadata",
                    str(meta_path),
                    "--exposure",
                    "group",
                    "--reference",
                    "ctrl",
                    "--covariates",
                    "sex",
                    "age",
                    "--output_dir",
                    str(out_dir),
                ]
            )
            self.assertEqual(res["tool"], "mg_multivariable")
            self.assertGreaterEqual(res["n_significant"], 1)
            self.assertEqual(res["top_hits"][0]["feature"], "Sp_pos")
            self.assertTrue((out_dir / "multivariable-clr.tsv").is_file())
            self.assertTrue((out_dir / "multivariable-volcano.png").is_file())

    def test_mg_mag_metabolism_ko_gapseq_sbml(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            ko_path = tmp_dir / "mag-ko-matrix.tsv"
            ko_path.write_text(
                "bin\tK02588\tK02586\tK02591\tK00249\tK01692\n"
                "bin_1\t1\t1\t1\t1\t0\n"
                "bin_2\t0\t0\t0\t1\t1\n",
                encoding="utf-8",
            )
            gp_path = tmp_dir / "bin_1-Pathways.tbl"
            gp_path.write_text(
                "ID\tName\tPrediction\tCompleteness\n"
                "PWY-1\tGlycolysis\ttrue\t90\n"
                "PWY-2\tMethanogenesis\tfalse\t10\n",
                encoding="utf-8",
            )
            sbml_path = tmp_dir / "bin_1.xml"
            sbml_path.write_text(
                "<?xml version='1.0' encoding='UTF-8'?>"
                "<sbml xmlns='http://www.sbml.org/sbml/level3/version1/core' xmlns:fbc='http://www.sbml.org/sbml/level3/version1/fbc/version2'>"
                "<model id='bin_1'>"
                "<listOfSpecies><species id='M_glc'/><species id='M_pyr'/></listOfSpecies>"
                "<listOfReactions><reaction id='R_EX_glc'/><reaction id='R_HEX1'/></listOfReactions>"
                "<fbc:listOfGeneProducts><fbc:geneProduct fbc:id='G_1'/></fbc:listOfGeneProducts>"
                "</model></sbml>",
                encoding="utf-8",
            )
            out_dir = tmp_dir / "out_met"
            res = _run_json(
                [
                    sys.executable,
                    str(MAG_SCRIPTS / "mag_metabolism.py"),
                    "--ko_matrix",
                    str(ko_path),
                    "--gapseq_pathways",
                    str(gp_path),
                    "--sbml",
                    str(sbml_path),
                    "--output_dir",
                    str(out_dir),
                ]
            )
            self.assertEqual(res["tool"], "mg_mag_metabolism")
            self.assertEqual(res["n_bins"], 2)
            self.assertEqual(res["n_sbml_models"], 1)
            self.assertTrue((out_dir / "mag-metabolism-summary.tsv").is_file())
            self.assertTrue((out_dir / "mag-module-completeness.tsv").is_file())
            self.assertTrue((out_dir / "mag-metabolism-heatmap.png").is_file())


if __name__ == "__main__":
    unittest.main()
