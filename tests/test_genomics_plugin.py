#!/usr/bin/env python3
"""End-to-end runtime tests for the genomics plugin (0.1.0):
- gen_inspect (MAF validation, VAF, context, clinical & CNV alignment)
- gen_summary (gene/sample summaries, Ti/Tv spectrum, dashboard & waterfall oncoplot)
- gen_tmb_ith (FoCR TMB with explicit capture_size_mb, MATH with CNV exclusion & group test)
- gen_interactions (pairwise co-occurrence/exclusivity, 10 Sanchez-Vega pathways, clinical contrast)
- gen_signatures (SBS96 catalog + NNLS refitting with cosine similarity + 6-class fallback)
- gen_cnv (CBS/CNVkit FGA, recurrent chromosome gain/loss, clinical group test)
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(
    REPO_ROOT,
    "resources",
    "plugins",
    "genomics",
    "skills",
    "genomics-workflow",
    "scripts",
)


def run_tool(script_name, args):
    cmd = [sys.executable, os.path.join(SCRIPTS_DIR, script_name)] + args
    env = {**os.environ, "MPLBACKEND": "Agg"}
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=180)
    if proc.returncode != 0:
        raise AssertionError(
            f"{script_name} failed (rc={proc.returncode})\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


class TestGenomicsPlugin(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.work = cls.tmpdir.name
        rng = np.random.default_rng(77)

        # 20 samples: 10 responder, 10 nonresponder
        samples = [f"S{i:02d}" for i in range(1, 21)]
        groups = ["responder"] * 10 + ["nonresponder"] * 10
        cls.clin_path = os.path.join(cls.work, "clinical.tsv")
        pd.DataFrame({"Tumor_Sample_Barcode": samples, "response": groups}).to_csv(
            cls.clin_path, sep="\t", index=False
        )

        # Build synthetic MAF with clear biological patterns:
        # - TP53 + PIK3CA co-occurring in responders
        # - KRAS mutually exclusive with EGFR
        # - Rich SBS1 (N[C>T]G) in responders vs SBS4 (C>A) in nonresponders
        maf_rows = []
        bases = ["A", "C", "G", "T"]
        for idx, s in enumerate(samples):
            is_resp = idx < 10
            # Driver hits
            if is_resp:
                for gene, vc in [("TP53", "Missense_Mutation"), ("PIK3CA", "Missense_Mutation"), ("EGFR", "Missense_Mutation")]:
                    if idx < 8:
                        maf_rows.append({
                            "Hugo_Symbol": gene, "Chromosome": "17", "Start_Position": 7577120 + idx,
                            "End_Position": 7577120 + idx, "Reference_Allele": "C", "Tumor_Seq_Allele2": "T",
                            "Variant_Classification": vc, "Variant_Type": "SNP",
                            "Tumor_Sample_Barcode": s, "t_vaf": round(float(rng.uniform(0.25, 0.55)), 4),
                            "Context": "ACG",
                        })
            else:
                if idx < 18:
                    maf_rows.append({
                        "Hugo_Symbol": "KRAS", "Chromosome": "12", "Start_Position": 25398284 + idx,
                        "End_Position": 25398284 + idx, "Reference_Allele": "C", "Tumor_Seq_Allele2": "A",
                        "Variant_Classification": "Missense_Mutation", "Variant_Type": "SNP",
                        "Tumor_Sample_Barcode": s, "t_vaf": round(float(rng.uniform(0.15, 0.50)), 4),
                        "Context": "CCA",
                    })
                    maf_rows.append({
                        "Hugo_Symbol": "APC", "Chromosome": "5", "Start_Position": 112175000 + idx,
                        "End_Position": 112175000 + idx, "Reference_Allele": "C", "Tumor_Seq_Allele2": "T",
                        "Variant_Classification": "Nonsense_Mutation", "Variant_Type": "SNP",
                        "Tumor_Sample_Barcode": s, "t_vaf": round(float(rng.uniform(0.10, 0.60)), 4),
                        "Context": "TCG",
                    })

            # Background SNVs (30 in responders -> higher TMB, 15 in nonresponders)
            n_bg = 30 if is_resp else 15
            for k in range(n_bg):
                b5 = bases[int(rng.integers(0, 4))]
                if is_resp and k % 2 == 0:
                    ref, alt, ctx = "C", "T", f"{b5}CG"
                else:
                    ref, alt, ctx = "C", "A", f"{b5}C{bases[int(rng.integers(0, 4))]}"
                vc = "Silent" if k % 5 == 0 else "Missense_Mutation"
                maf_rows.append({
                    "Hugo_Symbol": f"GENE{k:02d}", "Chromosome": str((k % 5) + 1),
                    "Start_Position": 100000 + k * 5000, "End_Position": 100000 + k * 5000,
                    "Reference_Allele": ref, "Tumor_Seq_Allele2": alt,
                    "Variant_Classification": vc, "Variant_Type": "SNP",
                    "Tumor_Sample_Barcode": s, "t_vaf": round(float(rng.uniform(0.08, 0.55)), 4),
                    "Context": ctx,
                })

        cls.maf_path = os.path.join(cls.work, "cohort.maf")
        pd.DataFrame(maf_rows).to_csv(cls.maf_path, sep="\t", index=False)

        # Synthetic 2-signature SBS96 reference matrix (SBS1_like: N[C>T]G, SBS4_like: N[C>A]N)
        sys.path.insert(0, SCRIPTS_DIR)
        from _common import all_sbs96_channels
        channels = all_sbs96_channels()
        sbs1 = np.array([1.0 if ch.endswith("[C>T]G") else 0.001 for ch in channels], dtype=float)
        sbs1 /= sbs1.sum()
        sbs4 = np.array([1.0 if "[C>A]" in ch else 0.001 for ch in channels], dtype=float)
        sbs4 /= sbs4.sum()
        cls.ref_sig_path = os.path.join(cls.work, "cosmic_ref.tsv")
        pd.DataFrame({"Type": channels, "SBS1": sbs1, "SBS4": sbs4}).to_csv(
            cls.ref_sig_path, sep="\t", index=False
        )

        # Synthetic CNV .seg file
        seg_rows = []
        for idx, s in enumerate(samples):
            is_resp = idx < 10
            seg_rows.append({"Sample": s, "Chromosome": "1", "Start": 1, "End": 10000000, "Segment_Mean": 0.02})
            seg_rows.append({"Sample": s, "Chromosome": "8", "Start": 1, "End": 8000000, "Segment_Mean": 0.65 if is_resp else 0.05})
            seg_rows.append({"Sample": s, "Chromosome": "17", "Start": 1, "End": 9000000, "Segment_Mean": -0.55 if not is_resp else -0.02})
        cls.seg_path = os.path.join(cls.work, "cohort.seg")
        pd.DataFrame(seg_rows).to_csv(cls.seg_path, sep="\t", index=False)

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def test_01_full_genomics_chain(self):
        # 1. gen_inspect
        insp = run_tool(
            "gen_inspect.py",
            ["--maf", self.maf_path, "--clinical", self.clin_path, "--seg", self.seg_path],
        )
        self.assertEqual(insp["n_samples"], 20)
        self.assertTrue(insp["vaf"]["present"])
        self.assertEqual(insp["context_column"], "Context")
        self.assertEqual(insp["clinical"]["n_matched_samples"], 20)
        self.assertEqual(insp["cnv_seg"]["n_matched_with_maf"], 20)

        # 2. gen_summary
        sum_dir = os.path.join(self.work, "summary")
        summ = run_tool("gen_summary.py", ["--maf", self.maf_path, "--top_genes", "15", "--output_dir", sum_dir])
        self.assertEqual(summ["n_samples"], 20)
        self.assertIsNotNone(summ["titv_ratio"])
        self.assertEqual(len(summ["figures"]), 2)

        # 3. gen_tmb_ith (with CNV exclusion + clinical group test)
        tmb_dir = os.path.join(self.work, "tmb_ith")
        tmb = run_tool(
            "gen_tmb_ith.py",
            [
                "--maf", self.maf_path,
                "--capture_size_mb", "2.0",
                "--min_vaf", "0.05",
                "--seg", self.seg_path,
                "--clinical", self.clin_path,
                "--group_col", "response",
                "--output_dir", tmb_dir,
            ],
        )
        self.assertTrue(tmb["cnv_filtered_for_math"])
        self.assertIsNotNone(tmb["tmb_median"])
        self.assertIsNotNone(tmb["math_median"])
        self.assertLess(tmb["group_comparison"]["tmb_per_mb"]["pvalue"], 0.01)

        # 4. gen_interactions (co-occurrence/exclusivity + 10 pathways + clinical contrast)
        int_dir = os.path.join(self.work, "interactions")
        inter = run_tool(
            "gen_interactions.py",
            [
                "--maf", self.maf_path,
                "--clinical", self.clin_path,
                "--group_col", "response",
                "--case", "responder",
                "--control", "nonresponder",
                "--top_genes", "15",
                "--output_dir", int_dir,
            ],
        )
        self.assertGreaterEqual(inter["n_significant_pairs"], 1)
        self.assertEqual(len(inter["pathways"]), 10)
        self.assertGreaterEqual(inter["cohort_comparison"]["n_significant"], 2)
        self.assertEqual(len(inter["figures"]), 2)

        # 5. gen_signatures (SBS96 + COSMIC NNLS refitting)
        sig_dir = os.path.join(self.work, "signatures")
        sig = run_tool(
            "gen_signatures.py",
            [
                "--maf", self.maf_path,
                "--ref_signatures", self.ref_sig_path,
                "--min_snvs", "10",
                "--output_dir", sig_dir,
            ],
        )
        self.assertEqual(sig["mode"], "sbs96")
        self.assertEqual(sig["refitting"]["n_samples_fitted"], 20)
        self.assertGreater(sig["refitting"]["median_cosine_similarity"], 0.75)
        self.assertEqual(len(sig["figures"]), 2)

        # 6. gen_cnv
        cnv_dir = os.path.join(self.work, "cnv")
        cnv = run_tool(
            "gen_cnv.py",
            [
                "--seg", self.seg_path,
                "--clinical", self.clin_path,
                "--group_col", "response",
                "--output_dir", cnv_dir,
            ],
        )
        self.assertEqual(cnv["n_samples"], 20)
        self.assertGreater(cnv["median_fga"], 0.2)
        self.assertEqual(len(cnv["figures"]), 1)


if __name__ == "__main__":
    unittest.main()
