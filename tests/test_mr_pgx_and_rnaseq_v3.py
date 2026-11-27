"""Contract, compilation, and end-to-end runtime tests for:
  - mendelian-randomization 0.1.0 (mr_gwas_qc, mr_harmonize, mr_estimate, mr_coloc)
  - pharmacogenomics 0.2.0 (pgx_dose_response, pgx_synergy, pgx_biomarker, pgx_hts_qc, pgx_predict_cohort)
  - bulk-rna-seq 0.4.0 upgrades (rna_coexpression, rna_tme_deconv, rna_regulon, rna_splicing)
"""
from __future__ import annotations

import ast
import json
import math
import py_compile
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PLUGINS_DIR = REPO_ROOT / "resources" / "plugins"
MR_DIR = PLUGINS_DIR / "mendelian-randomization"
PGX_DIR = PLUGINS_DIR / "pharmacogenomics"
RNA_DIR = PLUGINS_DIR / "bulk-rna-seq"

try:
    import matplotlib  # noqa: F401
    import numpy  # noqa: F401
    import pandas  # noqa: F401
    import scipy  # noqa: F401
    import sklearn  # noqa: F401

    HAS_SCIPY_STACK = True
except ImportError:
    HAS_SCIPY_STACK = False


def _run_json(cmd: list[str]) -> dict:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise AssertionError(
            f"Command failed ({proc.returncode}): {' '.join(cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


class TestPluginManifestsAndOfflineCompilation(unittest.TestCase):
    def test_mendelian_randomization_layout_and_compilation(self) -> None:
        manifest = (MR_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("id: mendelian-randomization", manifest)
        self.assertIn("version: 0.1.0", manifest)
        self.assertIn("toolPrefix: mr", manifest)
        self.assertTrue((MR_DIR / "icon.svg").is_file())
        self.assertTrue((MR_DIR / "README.md").is_file())
        self.assertTrue((MR_DIR / "CHANGELOG.md").is_file())
        self.assertTrue((MR_DIR / "agents" / "MendelianRandomization.md").is_file())

        locks_dir = MR_DIR / "environments" / "mr" / "locks"
        for plat in ("darwin-arm64.txt", "darwin-x64.txt", "linux-x64.txt"):
            self.assertTrue((locks_dir / plat).is_file(), f"Missing lock {plat}")

        skill_dir = MR_DIR / "skills" / "mr-workflow"
        for schema in ("gwas-qc-result.json", "harmonize-result.json", "estimate-result.json", "coloc-result.json"):
            obj = json.loads((skill_dir / "schemas" / schema).read_text(encoding="utf-8"))
            self.assertEqual(obj["type"], "object")

        for script in ("_common.py", "mr_gwas_qc.py", "mr_harmonize.py", "mr_estimate.py", "mr_coloc.py"):
            p = skill_dir / "scripts" / script
            self.assertTrue(p.is_file(), f"Missing script {script}")
            ast.parse(p.read_text(encoding="utf-8"), filename=str(p))
            with tempfile.TemporaryDirectory() as tmp:
                py_compile.compile(str(p), cfile=str(Path(tmp) / f"{script}.pyc"), doraise=True)

    def test_pharmacogenomics_layout_and_compilation(self) -> None:
        manifest = (PGX_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("id: pharmacogenomics", manifest)
        self.assertIn("version: 0.2.0", manifest)
        self.assertIn("toolPrefix: pgx", manifest)
        self.assertTrue((PGX_DIR / "icon.svg").is_file())
        self.assertTrue((PGX_DIR / "README.md").is_file())
        self.assertTrue((PGX_DIR / "CHANGELOG.md").is_file())
        self.assertTrue((PGX_DIR / "agents" / "Pharmacogenomics.md").is_file())

        locks_dir = PGX_DIR / "environments" / "pgx" / "locks"
        for plat in ("darwin-arm64.txt", "darwin-x64.txt", "linux-x64.txt"):
            self.assertTrue((locks_dir / plat).is_file(), f"Missing lock {plat}")

        skill_dir = PGX_DIR / "skills" / "pgx-workflow"
        skill_md = (skill_dir / "SKILL.md").read_text(encoding="utf-8")
        self.assertIn("- name: hts_qc", skill_md)
        self.assertIn("- name: predict_cohort", skill_md)
        for schema in (
            "dose-response-result.json",
            "synergy-result.json",
            "biomarker-result.json",
            "hts-qc-result.json",
            "predict-cohort-result.json",
        ):
            obj = json.loads((skill_dir / "schemas" / schema).read_text(encoding="utf-8"))
            self.assertEqual(obj["type"], "object")

        for script in (
            "_common.py",
            "pgx_dose_response.py",
            "pgx_synergy.py",
            "pgx_biomarker.py",
            "pgx_hts_qc.py",
            "pgx_predict_cohort.py",
        ):
            p = skill_dir / "scripts" / script
            self.assertTrue(p.is_file(), f"Missing script {script}")
            ast.parse(p.read_text(encoding="utf-8"), filename=str(p))
            with tempfile.TemporaryDirectory() as tmp:
                py_compile.compile(str(p), cfile=str(Path(tmp) / f"{script}.pyc"), doraise=True)

    def test_bulk_rnaseq_v3_layout_and_compilation(self) -> None:
        manifest = (RNA_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("version: 0.4.0", manifest)
        skill_dir = RNA_DIR / "skills" / "rnaseq-workflow"
        skill_md = (skill_dir / "SKILL.md").read_text(encoding="utf-8")
        self.assertIn("- name: coexpression", skill_md)
        self.assertIn("- name: tme_deconv", skill_md)
        self.assertIn("- name: regulon", skill_md)
        self.assertIn("- name: splicing", skill_md)

        for schema in (
            "coexpression-result.json",
            "tme-deconv-result.json",
            "regulon-result.json",
            "splicing-result.json",
        ):
            obj = json.loads((skill_dir / "schemas" / schema).read_text(encoding="utf-8"))
            self.assertEqual(obj["type"], "object")

        for script in (
            "_common.py",
            "rna_coexpression.py",
            "rna_tme_deconv.py",
            "rna_regulon.py",
            "rna_splicing.py",
        ):
            p = skill_dir / "scripts" / script
            self.assertTrue(p.is_file(), f"Missing script {script}")
            src = p.read_text(encoding="utf-8")
            self.assertNotIn("sklearn", src, f"bulk-rna-seq script {script} must not import sklearn")
            ast.parse(src, filename=str(p))
            with tempfile.TemporaryDirectory() as tmp:
                py_compile.compile(str(p), cfile=str(Path(tmp) / f"{script}.pyc"), doraise=True)


@unittest.skipUnless(HAS_SCIPY_STACK, "Requires numpy/pandas/scipy/scikit-learn/matplotlib")
class TestMendelianRandomizationRuntime(unittest.TestCase):
    def test_mr_pipeline_and_coloc_end_to_end(self) -> None:
        scripts = MR_DIR / "skills" / "mr-workflow" / "scripts"
        rng = random.Random(2026)
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            exp_rows = ["SNP\tCHR\tPOS\tEA\tNEA\tEAF\tBETA\tSE\tP\tN"]
            out_rows = ["SNP\tCHR\tPOS\tEA\tNEA\tEAF\tBETA\tSE\tP\tN"]
            true_causal = 0.45

            # 18 strong non-palindromic instruments + 1 flipped allele instrument + 1 ambiguous palindromic + 2 weak SNPs
            for i in range(22):
                snp = f"rs{10000 + i}"
                pos = 1_000_000 + i * 600_000  # well separated (>500 kb) for clumping
                if i < 18:
                    ea_e, nea_e = "A", "G"
                    eaf = round(0.18 + 0.015 * i, 3)
                    bx = 0.08 + 0.008 * i
                    sex = 0.011
                    ea_o, nea_o, eaf_o, sign_o = ea_e, nea_e, eaf, 1.0
                    if i == 5:
                        # Flipped effect/other alleles in outcome
                        ea_o, nea_o, eaf_o, sign_o = nea_e, ea_e, round(1.0 - eaf, 3), -1.0
                elif i == 18:
                    # Unambiguous palindromic A/T with EAF=0.22
                    ea_e, nea_e = "A", "T"
                    eaf = 0.22
                    bx, sex = 0.12, 0.012
                    ea_o, nea_o, eaf_o, sign_o = "T", "A", 0.78, -1.0
                elif i == 19:
                    # Ambiguous palindromic C/G with EAF=0.49 -> should be dropped in harmonization
                    ea_e, nea_e = "C", "G"
                    eaf = 0.49
                    bx, sex = 0.14, 0.012
                    ea_o, nea_o, eaf_o, sign_o = "C", "G", 0.49, 1.0
                else:
                    # Weak instrument (P ~ 0.2, F < 10)
                    ea_e, nea_e = "C", "T"
                    eaf = 0.30
                    bx, sex = 0.01, 0.02
                    ea_o, nea_o, eaf_o, sign_o = "C", "T", 0.30, 1.0

                from scipy import stats as sp_stats

                px = float(2.0 * sp_stats.norm.sf(abs(bx / sex)))
                by = sign_o * (true_causal * bx + rng.gauss(0.0, 0.004))
                sey = 0.015
                py = float(2.0 * sp_stats.norm.sf(abs(by / sey)))
                exp_rows.append(f"{snp}\t1\t{pos}\t{ea_e}\t{nea_e}\t{eaf}\t{bx:.5f}\t{sex:.5f}\t{px:.3e}\t50000")
                out_rows.append(f"{snp}\t1\t{pos}\t{ea_o}\t{nea_o}\t{eaf_o}\t{by:.5f}\t{sey:.5f}\t{py:.3e}\t60000")

            exp_path = tmp_dir / "exposure_gwas.tsv"
            out_path = tmp_dir / "outcome_gwas.tsv"
            exp_path.write_text("\n".join(exp_rows) + "\n", encoding="utf-8")
            out_path.write_text("\n".join(out_rows) + "\n", encoding="utf-8")

            qc_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "mr_gwas_qc.py"),
                    "--gwas",
                    str(exp_path),
                    "--clump_kb",
                    "500",
                    "--output_dir",
                    str(tmp_dir / "qc"),
                ]
            )
            self.assertEqual(qc_res["tool"], "mr_gwas_qc")
            self.assertGreaterEqual(qc_res["n_lead_loci"], 18)

            harm_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "mr_harmonize.py"),
                    "--exposure",
                    qc_res["outputs"]["lead_instruments"],
                    "--outcome",
                    str(out_path),
                    "--clump_kb",
                    "500",
                    "--output_dir",
                    str(tmp_dir / "harm"),
                ]
            )
            self.assertEqual(harm_res["tool"], "mr_harmonize")
            self.assertGreaterEqual(harm_res["n_harmonized_retained"], 18)
            self.assertGreaterEqual(harm_res["n_palindromic_dropped"], 1)

            est_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "mr_estimate.py"),
                    "--harmonized",
                    harm_res["outputs"]["harmonized_table"],
                    "--output_dir",
                    str(tmp_dir / "est"),
                ]
            )
            self.assertEqual(est_res["tool"], "mr_estimate")
            ivw_row = est_res["methods"][0]
            self.assertAlmostEqual(ivw_row["beta"], true_causal, delta=0.08)
            self.assertLess(ivw_row["pvalue"], 1e-6)
            self.assertTrue(est_res["steiger_directionality"]["correct_causal_direction"])
            self.assertTrue(Path(est_res["figures"][0] + ".phi-artifact.json").is_file())

            # Regional colocalization with shared causal variant at index 25
            reg1 = ["SNP\tCHR\tPOS\tEA\tNEA\tEAF\tBETA\tSE\tP\tN"]
            reg2 = ["SNP\tCHR\tPOS\tEA\tNEA\tEAF\tBETA\tSE\tP\tN"]
            for j in range(60):
                rs = f"rs9{j:04d}"
                p_loc = 5_000_000 + j * 2_000
                b1 = 0.28 if j == 25 else rng.gauss(0.0, 0.015)
                b2 = 0.24 if j == 25 else rng.gauss(0.0, 0.015)
                s1, s2 = 0.025, 0.025
                from scipy import stats as sp_stats

                p1 = float(2.0 * sp_stats.norm.sf(abs(b1 / s1)))
                p2 = float(2.0 * sp_stats.norm.sf(abs(b2 / s2)))
                reg1.append(f"{rs}\t2\t{p_loc}\tA\tG\t0.30\t{b1:.5f}\t{s1:.5f}\t{p1:.3e}\t40000")
                reg2.append(f"{rs}\t2\t{p_loc}\tA\tG\t0.30\t{b2:.5f}\t{s2:.5f}\t{p2:.3e}\t45000")

            r1_path = tmp_dir / "trait1_locus.tsv"
            r2_path = tmp_dir / "trait2_locus.tsv"
            r1_path.write_text("\n".join(reg1) + "\n", encoding="utf-8")
            r2_path.write_text("\n".join(reg2) + "\n", encoding="utf-8")

            coloc_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "mr_coloc.py"),
                    "--trait1",
                    str(r1_path),
                    "--trait2",
                    str(r2_path),
                    "--output_dir",
                    str(tmp_dir / "coloc"),
                ]
            )
            self.assertEqual(coloc_res["tool"], "mr_coloc")
            self.assertGreater(coloc_res["posteriors"]["PP.H4"], 0.85)
            self.assertEqual(coloc_res["top_snp"]["SNP"], "rs90025")


@unittest.skipUnless(HAS_SCIPY_STACK, "Requires numpy/pandas/scipy/scikit-learn/matplotlib")
class TestPharmacogenomicsRuntime(unittest.TestCase):
    def test_pgx_dose_response_synergy_and_biomarker(self) -> None:
        scripts = PGX_DIR / "skills" / "pgx-workflow" / "scripts"
        rng = random.Random(77)
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)

            # 1. Dose-response: sensitive line (IC50 ~ 0.25 uM) + resistant line (Einf = 0.72 -> beyond_tested_range)
            dr_rows = ["compound\tsample\tdose\tviability"]
            doses = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
            for d in doses:
                for rep in range(2):
                    y_sens = 100.0 * (0.05 + 0.95 / (1.0 + (d / 0.25) ** 1.5)) if d > 0 else 100.0
                    y_res = 100.0 * (0.72 + 0.28 / (1.0 + (d / 2.0) ** 1.2)) if d > 0 else 100.0
                    dr_rows.append(f"Olaparib\tCell_Sens\t{d}\t{y_sens + rng.gauss(0, 1.2):.2f}")
                    dr_rows.append(f"Olaparib\tCell_Res\t{d}\t{y_res + rng.gauss(0, 1.2):.2f}")
            dr_path = tmp_dir / "viability.tsv"
            dr_path.write_text("\n".join(dr_rows) + "\n", encoding="utf-8")

            dr_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "pgx_dose_response.py"),
                    "--input",
                    str(dr_path),
                    "--output_dir",
                    str(tmp_dir / "dr_out"),
                ]
            )
            self.assertEqual(dr_res["status"], "ok")
            self.assertEqual(dr_res["detected_viability_scale"], "percent")
            self.assertEqual(dr_res["n_curves_fitted"], 2)
            self.assertEqual(dr_res["n_beyond_tested_range"], 1)
            by_sample = {c["sample"]: c for c in dr_res["curves"]}
            self.assertFalse(by_sample["Cell_Sens"]["beyond_tested_range"])
            self.assertAlmostEqual(by_sample["Cell_Sens"]["ic50"], 0.25, delta=0.08)
            self.assertTrue(by_sample["Cell_Res"]["beyond_tested_range"])
            self.assertTrue(by_sample["Cell_Res"]["ic50_report"].startswith(">"))

            # 2. 2D Drug Synergy: strong excess inhibition at interior doses
            syn_rows = ["sample\tdrug_a\tdrug_b\tdose_a\tdose_b\tviability"]
            grid_d = [0.0, 0.1, 0.3, 1.0, 3.0]
            for da in grid_d:
                ea = 0.80 / (1.0 + (1.2 / da) ** 1.3) if da > 0 else 0.0
                for db in grid_d:
                    eb = 0.75 / (1.0 + (1.5 / db) ** 1.3) if db > 0 else 0.0
                    bliss = ea + eb - ea * eb
                    excess = 0.28 if (da >= 0.3 and db >= 0.3) else (0.15 if (da > 0 and db > 0) else 0.0)
                    e_tot = min(0.98, bliss + excess)
                    syn_rows.append(f"A375\tDabrafenib\tTrametinib\t{da}\t{db}\t{1.0 - e_tot:.4f}")
            syn_path = tmp_dir / "combo.tsv"
            syn_path.write_text("\n".join(syn_rows) + "\n", encoding="utf-8")

            syn_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "pgx_synergy.py"),
                    "--input",
                    str(syn_path),
                    "--output_dir",
                    str(tmp_dir / "syn_out"),
                ]
            )
            self.assertEqual(syn_res["status"], "ok")
            combo_info = syn_res["combinations"][0]
            self.assertEqual(combo_info["classification"], "synergistic")
            self.assertGreater(combo_info["mean_zip"], 10.0)
            self.assertGreater(combo_info["mean_bliss"], 10.0)

            # 3. Pharmacogenomic biomarker screening (MAF + LN_IC50 where lower = sensitive)
            samples = [f"CL_{i:02d}" for i in range(18)]
            maf_rows = ["Hugo_Symbol\tTumor_Sample_Barcode\tVariant_Classification"]
            drug_rows = ["sample\tPLX4720\tControlDrug"]
            for i, s in enumerate(samples):
                is_braf = i < 8
                if is_braf:
                    maf_rows.append(f"BRAF\t{s}\tMissense_Mutation")
                if i % 3 == 0:
                    maf_rows.append(f"TP53\t{s}\tNonsense_Mutation")
                ln_ic50 = (-2.2 + rng.gauss(0, 0.25)) if is_braf else (1.4 + rng.gauss(0, 0.3))
                ctrl = rng.gauss(0.5, 0.4)
                drug_rows.append(f"{s}\t{ln_ic50:.4f}\t{ctrl:.4f}")

            maf_path = tmp_dir / "somatic.maf"
            drug_path = tmp_dir / "gdsc_ln_ic50.tsv"
            maf_path.write_text("\n".join(maf_rows) + "\n", encoding="utf-8")
            drug_path.write_text("\n".join(drug_rows) + "\n", encoding="utf-8")

            bm_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "pgx_biomarker.py"),
                    "--drug_response",
                    str(drug_path),
                    "--features",
                    str(maf_path),
                    "--output_dir",
                    str(tmp_dir / "bm_out"),
                ]
            )
            self.assertEqual(bm_res["status"], "ok")
            self.assertEqual(bm_res["detected_feature_type"], "maf")
            top_hit = bm_res["top_biomarkers"][0]
            self.assertEqual(top_hit["feature"], "BRAF")
            self.assertEqual(top_hit["drug"], "PLX4720")
            self.assertEqual(top_hit["direction"], "sensitivity")
            self.assertGreater(top_hit["sensitivity_effect"], 2.0)
            self.assertLess(top_hit["padj"], 0.01)

    def test_pgx_hts_qc_and_predict_cohort(self) -> None:
        scripts = PGX_DIR / "skills" / "pgx-workflow" / "scripts"
        rng = random.Random(202602)
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)

            # 1. HTS 96-well plate with row/edge bias + 2 potent inhibitor hits
            plate_rows = ["plate\twell\twell_type\tsignal"]
            for r_idx in range(8):
                row_ch = chr(ord("A") + r_idx)
                row_drift = (r_idx - 3.5) * 180.0
                for c_idx in range(1, 13):
                    well = f"{row_ch}{c_idx:02d}"
                    if c_idx == 1:
                        wtype = "pos_ctrl"
                        sig = 800.0 + rng.gauss(0, 35.0)
                    elif c_idx == 12:
                        wtype = "neg_ctrl"
                        sig = 10000.0 + rng.gauss(0, 140.0)
                    else:
                        wtype = "sample"
                        sig = 10000.0 + row_drift + rng.gauss(0, 120.0)
                        if well in ("D05", "F08"):
                            sig = 1400.0 + row_drift
                    plate_rows.append(f"Plate_A\t{well}\t{wtype}\t{sig:.2f}")

            plate_path = tmp_dir / "hts_plates.tsv"
            plate_path.write_text("\n".join(plate_rows) + "\n", encoding="utf-8")

            hts_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "pgx_hts_qc.py"),
                    "--plates",
                    str(plate_path),
                    "--apply_bscore",
                    "--output_dir",
                    str(tmp_dir / "hts_out"),
                ]
            )
            self.assertEqual(hts_res["status"], "ok")
            self.assertEqual(hts_res["tool"], "pgx_hts_qc")
            self.assertTrue(hts_res["apply_bscore"])
            self.assertEqual(hts_res["n_passed_plates"], 1)
            self.assertGreater(hts_res["plate_summaries"][0]["z_prime"], 0.7)
            self.assertTrue(hts_res["plate_summaries"][0]["spatial_bias_detected"])
            self.assertGreaterEqual(hts_res["n_hits_total"], 2)
            hit_wells = {h["well"] for h in hts_res["top_hits"]}
            self.assertIn("D05", hit_wells)
            self.assertIn("F08", hit_wells)
            self.assertTrue(Path(hts_res["outputs"]["plate_qc_summary"] + ".phi-artifact.json").is_file())
            self.assertTrue(Path(hts_res["outputs"]["plate_heatmap_qc"] + ".phi-artifact.json").is_file())

            # 2. Cross-cohort drug sensitivity prediction (RidgeCV)
            n_lines, n_clin, n_genes = 36, 14, 40
            train_lines = [f"CL_{i:02d}" for i in range(n_lines)]
            clin_patients = [f"PT_{i:02d}" for i in range(n_clin)]
            lat_tr = [(-1.5 if i < 18 else 1.5) + rng.gauss(0, 0.2) for i in range(n_lines)]
            lat_cl = [(-1.4 if i < 7 else 1.4) + rng.gauss(0, 0.2) for i in range(n_clin)]

            tr_expr_lines = ["gene\t" + "\t".join(train_lines)]
            cl_expr_lines = ["gene\t" + "\t".join(clin_patients)]
            for g in range(n_genes):
                w = 1.6 if g < 15 else 0.0
                tr_vals = [7.0 + w * lat_tr[i] + rng.gauss(0, 0.4) for i in range(n_lines)]
                # Shifted baseline in clinical cohort to test cross-platform Z-score harmonization
                cl_vals = [10.5 + 1.4 * w * lat_cl[i] + rng.gauss(0, 0.5) for i in range(n_clin)]
                tr_expr_lines.append(f"G_{g:03d}\t" + "\t".join(f"{v:.4f}" for v in tr_vals))
                cl_expr_lines.append(f"G_{g:03d}\t" + "\t".join(f"{v:.4f}" for v in cl_vals))

            resp_lines = ["sample\tPredictableDrug\tNoiseDrug"]
            for i, s in enumerate(train_lines):
                y_pred = 1.2 * lat_tr[i] + rng.gauss(0, 0.25)
                y_noise = rng.gauss(0, 1.0)
                resp_lines.append(f"{s}\t{y_pred:.4f}\t{y_noise:.4f}")

            tr_expr_p = tmp_dir / "train_expr.tsv"
            cl_expr_p = tmp_dir / "clin_expr.tsv"
            resp_p = tmp_dir / "train_resp.tsv"
            tr_expr_p.write_text("\n".join(tr_expr_lines) + "\n", encoding="utf-8")
            cl_expr_p.write_text("\n".join(cl_expr_lines) + "\n", encoding="utf-8")
            resp_p.write_text("\n".join(resp_lines) + "\n", encoding="utf-8")

            pred_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "pgx_predict_cohort.py"),
                    "--train_expr",
                    str(tr_expr_p),
                    "--train_response",
                    str(resp_p),
                    "--clinical_expr",
                    str(cl_expr_p),
                    "--min_cv_r",
                    "0.30",
                    "--output_dir",
                    str(tmp_dir / "pred_out"),
                ]
            )
            self.assertEqual(pred_res["status"], "ok")
            self.assertEqual(pred_res["tool"], "pgx_predict_cohort")
            self.assertEqual(pred_res["n_drugs_modeled"], 2)
            self.assertGreaterEqual(pred_res["n_drugs_reliable"], 1)
            self.assertGreaterEqual(pred_res["n_drugs_low_predictability"], 1)
            by_drug = {m["drug"]: m for m in pred_res["drug_metrics"]}
            self.assertTrue(by_drug["PredictableDrug"]["reliable"])
            self.assertGreater(by_drug["PredictableDrug"]["cv_pearson_r"], 0.70)
            self.assertTrue(by_drug["NoiseDrug"]["low_predictability_warning"])
            self.assertTrue(Path(pred_res["outputs"]["predicted_response"] + ".phi-artifact.json").is_file())
            self.assertTrue(Path(pred_res["outputs"]["overview_plot"] + ".phi-artifact.json").is_file())


@unittest.skipUnless(HAS_SCIPY_STACK, "Requires numpy/pandas/scipy/scikit-learn/matplotlib")
class TestBulkRnaSeqV3Runtime(unittest.TestCase):
    def test_rna_coexpression_and_tme_deconv(self) -> None:
        scripts = RNA_DIR / "skills" / "rnaseq-workflow" / "scripts"
        rng = random.Random(314)
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)

            # 1. Co-expression cohort: 16 samples x 150 genes with 3 planted modules
            samples = [f"S{i:02d}" for i in range(16)]
            latent1 = [(-1.2 if i < 8 else 1.2) + rng.gauss(0, 0.25) for i in range(16)]
            latent2 = [math.sin(i * 0.5) + rng.gauss(0, 0.25) for i in range(16)]
            latent3 = [math.cos(i * 0.4) + rng.gauss(0, 0.25) for i in range(16)]

            expr_lines = ["gene\t" + "\t".join(samples)]
            for g in range(150):
                lat = latent1 if g < 50 else (latent2 if g < 100 else latent3)
                vals = [max(0.1, 6.0 + 1.8 * lat[s_i] + rng.gauss(0, 0.35)) for s_i in range(16)]
                expr_lines.append(f"GENE_{g:03d}\t" + "\t".join(f"{v:.4f}" for v in vals))

            expr_path = tmp_dir / "log_expr.tsv"
            expr_path.write_text("\n".join(expr_lines) + "\n", encoding="utf-8")

            meta_lines = ["sample\tcondition\tseverity_score"]
            for i, s in enumerate(samples):
                cond = "Responder" if i < 8 else "NonResponder"
                sev = 10.0 + 5.0 * latent1[i] + rng.gauss(0, 0.4)
                meta_lines.append(f"{s}\t{cond}\t{sev:.3f}")
            meta_path = tmp_dir / "meta.tsv"
            meta_path.write_text("\n".join(meta_lines) + "\n", encoding="utf-8")

            coex_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "rna_coexpression.py"),
                    "--input",
                    str(expr_path),
                    "--metadata",
                    str(meta_path),
                    "--min_module_size",
                    "15",
                    "--output_dir",
                    str(tmp_dir / "coex_out"),
                ]
            )
            self.assertEqual(coex_res["status"], "ok")
            self.assertGreaterEqual(coex_res["n_modules"], 2)
            self.assertTrue(len(coex_res["top_trait_associations"]) >= 2)
            self.assertLess(coex_res["top_trait_associations"][0]["padj"], 0.01)
            self.assertTrue(Path(coex_res["outputs"]["module_eigengenes"] + ".phi-artifact.json").is_file())

            # 2. TME Deconvolution: 10 samples (5 Hot CD8_T/Immune-high vs 5 Cold Stromal-high)
            if str(scripts) not in sys.path:
                sys.path.insert(0, str(scripts))
            from rna_tme_deconv import BUILTIN_IMMUNE_PANEL, IMMUNE_GENES, STROMAL_GENES  # type: ignore

            all_markers = sorted(
                set(STROMAL_GENES)
                | set(IMMUNE_GENES)
                | {g for prim, sec in BUILTIN_IMMUNE_PANEL.values() for g in (prim + sec)}
                | {f"HK_{k}" for k in range(60)}
            )
            tme_samples = [f"Hot_{i}" for i in range(5)] + [f"Cold_{i}" for i in range(5)]
            cd8_set = set(BUILTIN_IMMUNE_PANEL["CD8_T"][0] + BUILTIN_IMMUNE_PANEL["CD8_T"][1])
            m2_set = set(BUILTIN_IMMUNE_PANEL["Macrophage_M2"][0])
            stromal_set = set(STROMAL_GENES)
            immune_set = set(IMMUNE_GENES)

            tme_lines = ["gene\t" + "\t".join(tme_samples)]
            for g in all_markers:
                vals = []
                for idx in range(10):
                    is_hot = idx < 5
                    base = 25.0 + rng.uniform(0, 8.0)
                    if is_hot and (g in cd8_set or g in immune_set):
                        base += 280.0 if g in cd8_set else 120.0
                    elif (not is_hot) and (g in stromal_set or g in m2_set):
                        base += 220.0 if g in stromal_set else 140.0
                    vals.append(f"{base:.2f}")
                tme_lines.append(f"{g}\t" + "\t".join(vals))

            tme_expr = tmp_dir / "tme_tpm.tsv"
            tme_expr.write_text("\n".join(tme_lines) + "\n", encoding="utf-8")
            tme_meta = tmp_dir / "tme_meta.tsv"
            tme_meta.write_text(
                "sample\tphenotype\n"
                + "\n".join(f"{s}\t{'Hot' if i < 5 else 'Cold'}" for i, s in enumerate(tme_samples))
                + "\n",
                encoding="utf-8",
            )

            tme_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "rna_tme_deconv.py"),
                    "--input",
                    str(tme_expr),
                    "--metadata",
                    str(tme_meta),
                    "--group_col",
                    "phenotype",
                    "--output_dir",
                    str(tmp_dir / "tme_out"),
                ]
            )
            self.assertEqual(tme_res["status"], "ok")
            self.assertEqual(tme_res["n_samples"], 10)
            self.assertGreater(tme_res["mean_cell_fractions"]["CD8_T"], 0.15)
            self.assertTrue(len(tme_res["top_group_differences"]) >= 3)
            self.assertTrue(Path(tme_res["outputs"]["cell_fractions"] + ".phi-artifact.json").is_file())

    def test_rna_regulon_and_splicing(self) -> None:
        scripts = RNA_DIR / "skills" / "rnaseq-workflow" / "scripts"
        rng = random.Random(4242)
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)

            # 1. rna_regulon: built-in 14-pathway network (ULM matrix mode + MLM DE mode)
            if str(scripts) not in sys.path:
                sys.path.insert(0, str(scripts))
            from rna_regulon import BUILTIN_REGULON_NETWORK  # type: ignore

            all_genes = sorted({g for pairs in BUILTIN_REGULON_NETWORK.values() for g, _ in pairs} | {f"BG_{i}" for i in range(40)})
            hypoxia_map = dict(BUILTIN_REGULON_NETWORK["Hypoxia"])
            p53_map = dict(BUILTIN_REGULON_NETWORK["p53"])

            reg_samples = [f"H_{i}" for i in range(4)] + [f"N_{i}" for i in range(4)]
            mat_lines = ["gene\t" + "\t".join(reg_samples)]
            de_lines = ["gene\tlog2FoldChange\tstat\tpvalue\tpadj"]
            for g in all_genes:
                vals = []
                for s_i in range(8):
                    v = 6.0 + rng.gauss(0, 0.25)
                    if s_i < 4 and g in hypoxia_map:
                        v += 3.2 * hypoxia_map[g]
                    if s_i >= 4 and g in p53_map:
                        v += 2.6 * p53_map[g]
                    vals.append(f"{v:.3f}")
                mat_lines.append(f"{g}\t" + "\t".join(vals))

                stat_v = rng.gauss(0, 0.4)
                if g in hypoxia_map:
                    stat_v += 5.5 * hypoxia_map[g]
                de_lines.append(f"{g}\t{stat_v * 0.4:.3f}\t{stat_v:.3f}\t0.001\t0.01")

            mat_path = tmp_dir / "reg_expr.tsv"
            de_path = tmp_dir / "reg_de.tsv"
            mat_path.write_text("\n".join(mat_lines) + "\n", encoding="utf-8")
            de_path.write_text("\n".join(de_lines) + "\n", encoding="utf-8")

            ulm_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "rna_regulon.py"),
                    "--input",
                    str(mat_path),
                    "--method",
                    "ulm",
                    "--output_dir",
                    str(tmp_dir / "reg_ulm"),
                ]
            )
            self.assertEqual(ulm_res["status"], "ok")
            self.assertEqual(ulm_res["input_mode"], "matrix")
            self.assertEqual(ulm_res["n_regulons_tested"], 14)
            self.assertGreaterEqual(ulm_res["n_significant_fdr"], 2)
            self.assertTrue(Path(ulm_res["outputs"]["activities_table"] + ".phi-artifact.json").is_file())
            self.assertTrue(Path(ulm_res["outputs"]["heatmap_plot"] + ".phi-artifact.json").is_file())

            mlm_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "rna_regulon.py"),
                    "--input",
                    str(de_path),
                    "--method",
                    "mlm",
                    "--output_dir",
                    str(tmp_dir / "reg_mlm"),
                ]
            )
            self.assertEqual(mlm_res["status"], "ok")
            self.assertEqual(mlm_res["input_mode"], "de")
            self.assertEqual(mlm_res["top_regulons"][0]["source"], "Hypoxia")
            self.assertEqual(mlm_res["top_regulons"][0]["direction"], "activated")

            # 2. rna_splicing: rMATS SE + SUPPA2 dPSI combined
            rmats_rows = [
                "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\tPValue\tFDR\tIncLevel1\tIncLevel2\tIncLevelDifference",
                "1\tENSG1\tTP53\tchr17\t-\t1000\t1100\t500\t700\t1500\t1700\t1e-8\t1e-6\t0.85,0.82\t0.25,0.28\t0.57",
                "2\tENSG2\tFGFR2\tchr10\t-\t2000\t2150\t1600\t1800\t2500\t2700\t2e-6\t2e-4\t0.15,0.18\t0.70,0.72\t-0.54",
                "3\tENSG3\tGAPDH\tchr12\t+\t3000\t3090\t2600\t2800\t3400\t3600\t0.45\t0.72\t0.50,0.52\t0.48,0.50\t0.02",
            ]
            suppa_rows = [
                "event_id\tdPSI\tp-val",
                "TP53;RI:chr17:500:700-1500:1700:-\t-0.31\t1e-5",
                "MDM2;A5:chr12:100-200:250-350:+\t0.24\t4e-4",
            ]
            rmats_p = tmp_dir / "SE.MATS.JC.txt"
            suppa_p = tmp_dir / "events.dpsi"
            rmats_p.write_text("\n".join(rmats_rows) + "\n", encoding="utf-8")
            suppa_p.write_text("\n".join(suppa_rows) + "\n", encoding="utf-8")

            sp_res = _run_json(
                [
                    sys.executable,
                    str(scripts / "rna_splicing.py"),
                    "--splicing_inputs",
                    str(rmats_p),
                    str(suppa_p),
                    "--output_dir",
                    str(tmp_dir / "sp_out"),
                ]
            )
            self.assertEqual(sp_res["status"], "ok")
            self.assertEqual(sp_res["tool"], "rna_splicing")
            self.assertEqual(sp_res["n_input_files"], 2)
            self.assertEqual(sp_res["n_events_total"], 5)
            self.assertEqual(sp_res["n_events_significant"], 4)
            # TP53 SE exon_len = 100 (frameshift_candidate); FGFR2 SE exon_len = 150 (in_frame_modular); TP53 RI intron = 799 (frameshift_candidate)
            self.assertGreaterEqual(sp_res["n_frameshift_candidates"], 1)
            self.assertGreaterEqual(sp_res["n_in_frame_modular"], 1)
            self.assertGreaterEqual(sp_res["n_multi_hit_genes"], 1)
            self.assertEqual(sp_res["top_multi_hit_genes"][0]["gene"], "TP53")
            self.assertTrue(Path(sp_res["outputs"]["significant_events"] + ".phi-artifact.json").is_file())
            self.assertTrue(Path(sp_res["outputs"]["volcano_plot"] + ".phi-artifact.json").is_file())


if __name__ == "__main__":
    unittest.main()
