#!/usr/bin/env python3
"""Contract and end-to-end runtime tests for the genomics plugin (0.3.0):
- Contract checks: manifest 0.3.0, summary <= 260, 9 tools/schemas, py_compile, Genomics.md
- gen_inspect (MAF validation, VAF, context, clinical & CNV alignment)
- gen_summary (gene/sample summaries, Ti/Tv spectrum, dashboard & waterfall oncoplot)
- gen_tmb_ith (FoCR TMB with explicit capture_size_mb, MATH with CNV exclusion & group test)
- gen_interactions (pairwise co-occurrence/exclusivity, 10 Sanchez-Vega pathways, clinical contrast)
- gen_signatures (SBS96 catalog + NNLS refitting with cosine similarity + 6-class fallback)
- gen_cnv (CBS/CNVkit FGA, recurrent chromosome gain/loss, clinical group test)
- gen_hrd_msi (scarHRD LOH/TAI/LST + MSIsensor-pro & MAF indel surrogate stratification)
- gen_clonality (McGranahan/PyClone CCF multiplicity MLE, 95% CI, drivers & 1/f neutral tail)
"""

import importlib.util
import json
import os
import py_compile
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HAS_SCIENTIFIC_STACK = all(
    importlib.util.find_spec(mod) is not None for mod in ("numpy", "pandas", "scipy", "matplotlib")
)
if HAS_SCIENTIFIC_STACK:
    import numpy as np
    import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLUGIN_DIR = Path(REPO_ROOT) / "resources" / "plugins" / "genomics"
SKILL_DIR = PLUGIN_DIR / "skills" / "genomics-workflow"
SCRIPTS_DIR = str(SKILL_DIR / "scripts")


def run_tool(script_name, args):
    cmd = [sys.executable, os.path.join(SCRIPTS_DIR, script_name)] + args
    env = {**os.environ, "MPLBACKEND": "Agg"}
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=180)
    if proc.returncode != 0:
        raise AssertionError(
            f"{script_name} failed (rc={proc.returncode})\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


class TestGenomicsContract(unittest.TestCase):
    def test_manifest_docs_skills_and_schemas_v030(self):
        manifest = (PLUGIN_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertRegex(manifest, re.compile(r"^version:\s*0\.3\.0\s*$", re.M))
        m_sum = re.search(r"^summary:\s*>-\n((?:\s+.*\n)+)", manifest, re.M)
        self.assertIsNotNone(m_sum)
        summary_text = " ".join(line.strip() for line in m_sum.group(1).splitlines())
        self.assertLessEqual(len(summary_text), 260)

        readme = (PLUGIN_DIR / "README.md").read_text(encoding="utf-8")
        self.assertIn("0.3.0", readme)
        self.assertIn("9 tools", readme)
        self.assertIn("gen_hrd_msi", readme)
        self.assertIn("gen_clonality", readme)

        changelog = (PLUGIN_DIR / "CHANGELOG.md").read_text(encoding="utf-8")
        self.assertIn("## 0.3.0 (2026-11-26)", changelog)

        agent_md = (PLUGIN_DIR / "agents" / "Genomics.md").read_text(encoding="utf-8")
        for kw in (
            "gen_hrd_msi",
            "gen_clonality",
            "hrd",
            "loh",
            "lst",
            "tai",
            "msi",
            "clonality",
            "ccf",
            "subclonal",
            "同源重组缺陷",
            "微卫星不稳定",
            "克隆演化",
        ):
            self.assertIn(kw, agent_md)

        skill_md = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
        fm = skill_md.split("\n---\n", 1)[0]
        self.assertNotRegex(fm, re.compile(r"^\s+toolPrefix:", re.M))
        tool_names = re.findall(r"^    - name:\s+(\S+)", fm, re.M)
        output_schemas = re.findall(r"^      output:\s+(\S+)", fm, re.M)
        self.assertEqual(len(tool_names), 9)
        self.assertEqual(len(output_schemas), 9)
        self.assertIn("hrd_msi", tool_names)
        self.assertIn("clonality", tool_names)
        for rel in output_schemas:
            schema_path = SKILL_DIR / rel
            self.assertTrue(schema_path.is_file(), f"missing schema {rel}")
            json.loads(schema_path.read_text(encoding="utf-8"))

        playbook = (SKILL_DIR / "references" / "analysis-playbook.md").read_text(encoding="utf-8")
        self.assertIn("gen_hrd_msi", playbook)
        self.assertIn("gen_clonality", playbook)

        for py_file in sorted((SKILL_DIR / "scripts").glob("*.py")):
            py_compile.compile(str(py_file), doraise=True)


@unittest.skipUnless(HAS_SCIENTIFIC_STACK, "requires plugin:genomics scientific stack (numpy/pandas/scipy/matplotlib)")
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
        self.assertEqual(len(summ["artifacts"]), 5)
        for art in summ["artifacts"]:
            self.assertTrue(os.path.exists(art + ".phi-artifact.json"))

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

    def test_02_hrd_msi_and_clonality(self):
        # 7. gen_hrd_msi: exact LOH / TAI / LST verification on hg38
        ascat_path = os.path.join(self.work, "ascat_segments.tsv")
        ascat_rows = [
            # HRD_SAMPLE:
            # chr1 p-arm LOH+TAI (25.5 Mb > 15 Mb & >= 11 Mb near telomere, not whole-arm): +1 LOH, +1 TAI
            {"sample": "HRD_SAMPLE", "chrom": "chr1", "start": 500000, "end": 26000000, "nMajor": 2, "nMinor": 0},
            # chr1 p-arm adjacent balanced gain (16 Mb >= 10 Mb): +1 LST break with prev
            {"sample": "HRD_SAMPLE", "chrom": "chr1", "start": 26000001, "end": 42000000, "nMajor": 2, "nMinor": 1},
            # chr1 p-arm tiny < 3 Mb noise (1.2 Mb): smoothed out before LST
            {"sample": "HRD_SAMPLE", "chrom": "chr1", "start": 42000001, "end": 43200000, "nMajor": 1, "nMinor": 1},
            # chr1 p-arm gain (16.8 Mb >= 10 Mb): +1 LST break with 26M..42M after smoothing 1.2 Mb
            {"sample": "HRD_SAMPLE", "chrom": "chr1", "start": 43200001, "end": 60000000, "nMajor": 3, "nMinor": 1},
            # chr1 q-arm telomeric imbalance (28.5 Mb >= 11 Mb, within 2 Mb of chr1 end): +1 TAI
            {"sample": "HRD_SAMPLE", "chrom": "chr1", "start": 220000000, "end": 248500000, "nMajor": 3, "nMinor": 1},
            # chr2 whole-chromosome LOH: must be excluded from LOH, TAI, LST (0, 0, 0)
            {"sample": "HRD_SAMPLE", "chrom": "chr2", "start": 1, "end": 242000000, "nMajor": 2, "nMinor": 0},
            # chr3 whole p-arm imbalance (reaches centromere): must be excluded from TAI (0 TAI)
            {"sample": "HRD_SAMPLE", "chrom": "chr3", "start": 1, "end": 87000000, "nMajor": 2, "nMinor": 1},
            # chr3 interstitial q-arm LOH (25 Mb > 15 Mb, away from telomere): +1 LOH
            {"sample": "HRD_SAMPLE", "chrom": "chr3", "start": 120000000, "end": 145000000, "nMajor": 2, "nMinor": 0},
            # MSI_SAMPLE & HRP_MSS_SAMPLE: balanced diploid
            {"sample": "MSI_SAMPLE", "chrom": "chr1", "start": 1, "end": 248000000, "nMajor": 1, "nMinor": 1},
            {"sample": "HRP_MSS_SAMPLE", "chrom": "chr1", "start": 1, "end": 248000000, "nMajor": 1, "nMinor": 1},
        ]
        pd.DataFrame(ascat_rows).to_csv(ascat_path, sep="\t", index=False)

        msi_path = os.path.join(self.work, "msisensor.tsv")
        pd.DataFrame([
            {"sample": "HRD_SAMPLE", "Total_Number_of_Sites": 1000, "Number_of_Somatic_Sites": 8, "%": 0.8},
            {"sample": "MSI_SAMPLE", "Total_Number_of_Sites": 1000, "Number_of_Somatic_Sites": 185, "%": 18.5},
            {"sample": "HRP_MSS_SAMPLE", "Total_Number_of_Sites": 1000, "Number_of_Somatic_Sites": 4, "%": 0.4},
        ]).to_csv(msi_path, sep="\t", index=False)

        hrd_dir = os.path.join(self.work, "hrd_msi")
        hrd_res = run_tool(
            "gen_hrd_msi.py",
            [
                "--segments", ascat_path,
                "--msi_table", msi_path,
                "--genome", "hg38",
                "--hrd_cutoff", "5",
                "--msi_pct_cutoff", "3.5",
                "--output_dir", hrd_dir,
            ],
        )
        self.assertEqual(hrd_res["n_samples"], 3)
        by_sample = {r["sample"]: r for r in hrd_res["samples"]}
        self.assertEqual(by_sample["HRD_SAMPLE"]["loh_score"], 2)
        self.assertEqual(by_sample["HRD_SAMPLE"]["tai_score"], 2)
        self.assertEqual(by_sample["HRD_SAMPLE"]["lst_score"], 2)
        self.assertEqual(by_sample["HRD_SAMPLE"]["hrd_sum"], 6)
        self.assertTrue(by_sample["HRD_SAMPLE"]["hrd_positive"])
        self.assertEqual(by_sample["HRD_SAMPLE"]["biomarker_stratum"], "HRD_pos_MSS")
        self.assertEqual(by_sample["MSI_SAMPLE"]["msi_status"], "MSI-H")
        self.assertEqual(by_sample["MSI_SAMPLE"]["biomarker_stratum"], "MSI_H")
        self.assertEqual(by_sample["HRP_MSS_SAMPLE"]["biomarker_stratum"], "HR_proficient_MSS")
        for art in hrd_res["artifacts"]:
            self.assertTrue(os.path.exists(art + ".phi-artifact.json"))

        # 8. gen_clonality: purity=0.60, CN=3 -> clonal CCF~1.0 vs subclonal CCF~0.35
        pur_path = os.path.join(self.work, "purity.tsv")
        pd.DataFrame([{"sample": "CLONAL_TUMOR", "purity": 0.60, "ploidy": 3.0}]).to_csv(
            pur_path, sep="\t", index=False
        )
        cn_seg_path = os.path.join(self.work, "clonality_segs.tsv")
        pd.DataFrame([
            {"sample": "CLONAL_TUMOR", "chrom": "17", "start": 1, "end": 80000000, "total_cn": 3},
            {"sample": "CLONAL_TUMOR", "chrom": "3", "start": 1, "end": 190000000, "total_cn": 3},
            {"sample": "CLONAL_TUMOR", "chrom": "12", "start": 1, "end": 130000000, "total_cn": 3},
        ]).to_csv(cn_seg_path, sep="\t", index=False)

        # Under purity=0.60, CN=3: f_exp(1) = 0.60 / (0.60*3 + 0.40*2) = 0.23077, f_exp(2) = 0.46154
        clon_rows = [
            # Clonal multiplicity=1 (VAF=0.23 -> CCF ~ 0.9967)
            {"Hugo_Symbol": "TP53", "Tumor_Sample_Barcode": "CLONAL_TUMOR", "Chromosome": "17", "Start_Position": 7577120, "Variant_Classification": "Missense_Mutation", "t_alt_count": 23, "t_ref_count": 77},
            # Clonal multiplicity=2 (VAF=0.46 -> CCF ~ 0.9967)
            {"Hugo_Symbol": "PIK3CA", "Tumor_Sample_Barcode": "CLONAL_TUMOR", "Chromosome": "3", "Start_Position": 178936091, "Variant_Classification": "Missense_Mutation", "t_alt_count": 46, "t_ref_count": 54},
            # Subclonal CCF ~ 0.35 (VAF=0.08 -> CCF = 0.08 * 2.60 / 0.60 = 0.3467)
            {"Hugo_Symbol": "KRAS", "Tumor_Sample_Barcode": "CLONAL_TUMOR", "Chromosome": "12", "Start_Position": 25398284, "Variant_Classification": "Missense_Mutation", "t_alt_count": 8, "t_ref_count": 92},
        ]
        # Add neutral-tail variants in [0.12, 0.24] following M(f) ~ 1/f
        for idx, alt_c in enumerate([12, 13, 14, 15, 16, 18, 20, 22, 24], start=1):
            clon_rows.append({
                "Hugo_Symbol": f"PASSENGER{idx}",
                "Tumor_Sample_Barcode": "CLONAL_TUMOR",
                "Chromosome": "12",
                "Start_Position": 50000000 + idx * 10000,
                "Variant_Classification": "Missense_Mutation",
                "t_alt_count": alt_c,
                "t_ref_count": 100 - alt_c,
            })

        clon_maf = os.path.join(self.work, "clonality.maf")
        pd.DataFrame(clon_rows).to_csv(clon_maf, sep="\t", index=False)

        clon_dir = os.path.join(self.work, "clonality")
        clon_res = run_tool(
            "gen_clonality.py",
            [
                "--maf", clon_maf,
                "--purity_table", pur_path,
                "--segments", cn_seg_path,
                "--output_dir", clon_dir,
            ],
        )
        self.assertEqual(clon_res["n_samples"], 1)
        var_df = pd.read_csv(clon_res["variants_table"], sep="\t")
        by_gene = {r["Hugo_Symbol"]: r for _, r in var_df.iterrows()}

        self.assertEqual(by_gene["TP53"]["multiplicity"], 1)
        self.assertAlmostEqual(by_gene["TP53"]["ccf"], 0.9967, places=2)
        self.assertEqual(by_gene["TP53"]["clonality"], "clonal")
        self.assertTrue(bool(by_gene["TP53"]["is_driver"]))

        self.assertEqual(by_gene["PIK3CA"]["multiplicity"], 2)
        self.assertAlmostEqual(by_gene["PIK3CA"]["ccf"], 0.9967, places=2)
        self.assertEqual(by_gene["PIK3CA"]["clonality"], "clonal")

        self.assertEqual(by_gene["KRAS"]["multiplicity"], 1)
        self.assertAlmostEqual(by_gene["KRAS"]["ccf"], 0.3467, places=2)
        self.assertEqual(by_gene["KRAS"]["clonality"], "subclonal")
        self.assertIsNotNone(clon_res["median_r2_neutral"])
        self.assertGreater(clon_res["median_r2_neutral"], 0.80)
        for art in clon_res["artifacts"]:
            self.assertTrue(os.path.exists(art + ".phi-artifact.json"))


if __name__ == "__main__":
    unittest.main()
