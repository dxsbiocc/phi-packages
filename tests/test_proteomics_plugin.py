#!/usr/bin/env python3
"""End-to-end runtime tests for the proteomics plugin (0.1.0):
- prot_inspect (MaxQuant-style LFQ table, contaminants/reverse, raw scale, design & batch collinearity)
- prot_qc_impute (contaminant removal, log2, median norm, MinProb & kNN imputation, PCA/QC plots)
- prot_de (empirical-Bayes moderated t-test + volcano/MA)
- prot_ptm (total-protein abundance normalization + KSEA kinase activity)
- prot_enrichment (detected-proteome ORA + deterministic pre-ranked GSEA)
- prot_network (Spearman co-abundance network, modules, hubs, heatmap)
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
    "proteomics",
    "skills",
    "proteomics-workflow",
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


class TestProteomicsPlugin(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.work = cls.tmpdir.name
        rng = np.random.default_rng(123)

        # 8 samples: 4 control, 4 treated across 2 balanced batches
        samples = [f"CTRL_{i}" for i in range(1, 5)] + [f"TRT_{i}" for i in range(1, 5)]
        cls.meta_path = os.path.join(cls.work, "metadata.tsv")
        meta = pd.DataFrame(
            {
                "sample": samples,
                "condition": ["control"] * 4 + ["treated"] * 4,
                "batch": ["B1", "B1", "B2", "B2", "B1", "B1", "B2", "B2"],
            }
        )
        meta.to_csv(cls.meta_path, sep="\t", index=False)

        # Synthetic MaxQuant-style raw LFQ intensity table (80 real proteins + 3 contaminants/decoys)
        n_real = 80
        genes = [f"PROT{i:03d}" for i in range(n_real)]
        # Realistic proteomics: wide dynamic range across proteins, tight replicate SD (~0.25)
        log2_mat = rng.normal(loc=24.0, scale=2.0, size=(n_real, 1)) + rng.normal(
            loc=0.0, scale=0.25, size=(n_real, 8)
        )
        # Strong upregulated module (PROT000..PROT011) and downregulated module (PROT012..PROT023)
        log2_mat[:12, 4:] += 2.6 + rng.normal(0, 0.15, size=(12, 4))
        log2_mat[12:24, 4:] -= 2.4 + rng.normal(0, 0.15, size=(12, 4))
        # Inject realistic MNAR missingness (zeros in raw LFQ)
        raw_lfq = np.power(2.0, log2_mat)
        raw_lfq[0, 0] = 0.0
        raw_lfq[1, 1] = 0.0
        raw_lfq[25, 0] = 0.0
        raw_lfq[25, 5] = 0.0

        rows = []
        for idx, g in enumerate(genes):
            row = {
                "Majority protein IDs": g,
                "Genes": g,
                "Potential contaminant": "",
                "Reverse": "",
                "Only identified by site": "",
            }
            for s_idx, s in enumerate(samples):
                row[f"LFQ intensity {s}"] = round(float(raw_lfq[idx, s_idx]), 2)
            rows.append(row)
        # Add 3 contaminant / decoy rows
        for cid, pc, rev in [("CON__KRT1", "+", ""), ("REV__P0001", "", "+"), ("SITE_ONLY", "", "")]:
            row = {
                "Majority protein IDs": cid,
                "Genes": cid,
                "Potential contaminant": pc,
                "Reverse": rev,
                "Only identified by site": "+" if cid == "SITE_ONLY" else "",
            }
            for s in samples:
                row[f"LFQ intensity {s}"] = 1e7
            rows.append(row)

        cls.mq_path = os.path.join(cls.work, "proteinGroups.txt")
        pd.DataFrame(rows).to_csv(cls.mq_path, sep="\t", index=False)

        # GMT gene-set file
        cls.gmt_path = os.path.join(cls.work, "pathways.gmt")
        with open(cls.gmt_path, "w", encoding="utf-8") as handle:
            handle.write("UP_SIGNALING\tNA\t" + "\t".join(genes[:12]) + "\n")
            handle.write("DOWN_METABOLISM\tNA\t" + "\t".join(genes[12:24]) + "\n")
            handle.write("BACKGROUND_COMPLEX\tNA\t" + "\t".join(genes[30:45]) + "\n")

        # Synthetic PTM/phosphosite matrix (log2) + kinase-substrate table
        ptm_rows = []
        kinase_rows = []
        for i in range(40):
            parent = genes[i]
            site = f"{parent}_S{100 + i}"
            base = log2_mat[i].copy()
            # Sites 0..7 have genuine stoichiometry increase beyond protein abundance
            if i < 8:
                base[4:] += 2.2
                kinase_rows.append({"kinase": "MAPK1", "site_id": site})
            elif i < 16:
                base[4:] -= 2.0
                kinase_rows.append({"kinase": "CDK2", "site_id": site})
            elif i < 24:
                kinase_rows.append({"kinase": "AKT1", "site_id": site})
            entry = {"site_id": site, "protein_id": parent}
            for s_idx, s in enumerate(samples):
                entry[s] = round(float(base[s_idx] + rng.normal(0, 0.15)), 4)
            ptm_rows.append(entry)

        cls.ptm_path = os.path.join(cls.work, "phospho_sites.tsv")
        pd.DataFrame(ptm_rows).to_csv(cls.ptm_path, sep="\t", index=False)
        cls.kinase_path = os.path.join(cls.work, "kinase_substrates.tsv")
        pd.DataFrame(kinase_rows).to_csv(cls.kinase_path, sep="\t", index=False)

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def test_01_inspect_qc_de_enrichment_network_and_ptm(self):
        # 1. prot_inspect
        insp = run_tool(
            "prot_inspect.py",
            [
                "--input", self.mq_path,
                "--metadata", self.meta_path,
                "--group_col", "condition",
                "--batch_col", "batch",
            ],
        )
        self.assertEqual(insp["detected_format"], "maxquant")
        self.assertEqual(insp["contaminant_or_decoy_rows"], 3)
        self.assertEqual(insp["scale"]["scale"], "raw_intensity")
        self.assertEqual(insp["design"]["matched_samples"], 8)
        self.assertFalse(insp["design"]["batch_collinearity"]["confounded"])

        # 2. prot_qc_impute (minprob)
        qc_dir = os.path.join(self.work, "qc")
        qc = run_tool(
            "prot_qc_impute.py",
            [
                "--input", self.mq_path,
                "--metadata", self.meta_path,
                "--normalization", "median",
                "--imputation", "minprob",
                "--output_dir", qc_dir,
            ],
        )
        self.assertEqual(qc["n_contaminants_removed"], 3)
        self.assertEqual(qc["n_features_kept"], 80)
        self.assertTrue(qc["log2_transformed"])
        self.assertEqual(len(qc["figures"]), 3)
        proc_mat = qc["processed_matrix"]
        self.assertTrue(os.path.exists(proc_mat))

        # 3. prot_de (empirical-Bayes moderated t-test)
        de_dir = os.path.join(self.work, "de")
        de = run_tool(
            "prot_de.py",
            [
                "--input", proc_mat,
                "--metadata", self.meta_path,
                "--case", "treated",
                "--control", "control",
                "--method", "moderated_t",
                "--output_dir", de_dir,
            ],
        )
        self.assertGreaterEqual(de["n_up"], 10)
        self.assertGreaterEqual(de["n_down"], 10)
        self.assertTrue(de["empirical_bayes"]["shrunk"])
        self.assertEqual(len(de["figures"]), 2)
        de_table = de["de_table"]

        # 4. prot_ptm (with total-protein normalization + KSEA)
        ptm_dir = os.path.join(self.work, "ptm")
        ptm = run_tool(
            "prot_ptm.py",
            [
                "--ptm_matrix", self.ptm_path,
                "--protein_matrix", proc_mat,
                "--metadata", self.meta_path,
                "--case", "treated",
                "--control", "control",
                "--kinase_table", self.kinase_path,
                "--min_substrates", "3",
                "--output_dir", ptm_dir,
            ],
        )
        self.assertTrue(ptm["protein_normalized"])
        self.assertEqual(ptm["n_sites_protein_normalized"], 40)
        self.assertGreaterEqual(ptm["n_up"], 6)
        self.assertEqual(ptm["ksea"]["n_kinases_tested"], 3)
        self.assertGreaterEqual(ptm["ksea"]["n_significant"], 1)

        # 5. prot_enrichment (ORA + GSEA)
        enr_dir = os.path.join(self.work, "enrichment")
        enr = run_tool(
            "prot_enrichment.py",
            [
                "--de_table", de_table,
                "--gmt", self.gmt_path,
                "--mode", "both",
                "--permutations", "200",
                "--output_dir", enr_dir,
            ],
        )
        self.assertEqual(enr["ora"]["universe_size"], 80)
        self.assertGreaterEqual(enr["ora"]["n_significant"], 1)
        self.assertGreaterEqual(enr["gsea"]["n_significant"], 1)
        self.assertEqual(len(enr["figures"]), 1)

        # 6. prot_network (co-abundance modules & hubs)
        net_dir = os.path.join(self.work, "network")
        net = run_tool(
            "prot_network.py",
            [
                "--input", proc_mat,
                "--de_table", de_table,
                "--top_n", "40",
                "--min_abs_corr", "0.75",
                "--output_dir", net_dir,
            ],
        )
        self.assertGreater(net["n_edges"], 10)
        self.assertGreaterEqual(net["n_modules_ge2"], 1)
        self.assertEqual(len(net["figures"]), 1)


if __name__ == "__main__":
    unittest.main()
