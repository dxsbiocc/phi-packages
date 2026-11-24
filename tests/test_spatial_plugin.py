#!/usr/bin/env python3
"""End-to-end runtime tests for spatial-transcriptomics 0.2.0 upgrades:
- st_py_spatial (SVG, nhood enrichment, co-occurrence, interaction_matrix,
  centrality_scores, Ripley's L, cross-sample spatial_diff)
- st_py_niche (composition + expression modes, cross-condition niche dynamics)
- st_py_communication (global + stratify_by niche-specific ligrec)
- st_py_trajectory (PAGA + spatial DPT + score_genes gene programs)
- st_py_notebook (rendering + self_check execution for niche & trajectory)
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(
    REPO_ROOT,
    "resources",
    "plugins",
    "spatial-transcriptomics",
    "skills",
    "spatial-workflow",
    "scripts",
)


def run_tool(script_name, args):
    cmd = [sys.executable, os.path.join(SCRIPTS_DIR, script_name)] + args
    env = {**os.environ, "MPLBACKEND": "Agg"}
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)
    if proc.returncode != 0:
        raise AssertionError(
            f"{script_name} failed (rc={proc.returncode})\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


def make_synthetic_spatial_h5ad(path):
    """Create a 120-spot, 60-gene clustered spatial AnnData with 3 domains and 2 conditions."""
    rng = np.random.default_rng(42)
    nx, ny = 12, 10
    n_obs = nx * ny
    gx, gy = np.meshgrid(np.arange(nx, dtype=float), np.arange(ny, dtype=float))
    coords = np.column_stack([gx.ravel() * 10.0, gy.ravel() * 10.0])

    # Three spatial regions along the x-axis.
    region = np.where(gx.ravel() < 4, "0", np.where(gx.ravel() < 8, "1", "2"))
    # Two conditions along the y-axis with slight region shift.
    condition = np.where(gy.ravel() < 5, "control", "treated")

    lr_genes = [
        "TGFB1", "TGFBR1", "CCL2", "CCR2", "VEGFA", "KDR",
        "SPP1", "CD44", "FN1", "ITGB1", "IL6", "IL6R",
    ]
    extra_genes = [f"GENE{i:02d}" for i in range(48)]
    genes = lr_genes + extra_genes
    n_vars = len(genes)

    counts = rng.poisson(lam=2.0, size=(n_obs, n_vars)).astype(np.float32)
    # Inject domain-specific expression + a smooth continuous spatial gradient along x
    # so the expression manifold connects domains 0 -> 1 -> 2 for PAGA/DPT.
    r0 = region == "0"
    r1 = region == "1"
    r2 = region == "2"
    x_norm = gx.ravel() / float(nx - 1)
    counts[r0, 0] += rng.poisson(25, size=r0.sum())   # TGFB1 high in domain 0
    counts[r1, 1] += rng.poisson(25, size=r1.sum())   # TGFBR1 high in domain 1
    counts[r1, 2] += rng.poisson(22, size=r1.sum())   # CCL2 high in domain 1
    counts[r2, 3] += rng.poisson(22, size=r2.sum())   # CCR2 high in domain 2
    counts[r0, 4] += rng.poisson(18, size=r0.sum())   # VEGFA
    counts[r2, 5] += rng.poisson(18, size=r2.sum())   # KDR
    for j in range(12, 24):
        counts[:, j] += rng.poisson(np.clip(12.0 * (1.0 - x_norm), 0.1, 20.0))
    for j in range(24, 36):
        counts[:, j] += rng.poisson(np.clip(12.0 * x_norm, 0.1, 20.0))

    adata = ad.AnnData(
        X=counts.copy(),
        obs=pd.DataFrame(
            {
                "leiden": pd.Categorical(region),
                "condition": pd.Categorical(condition),
                "sample": pd.Categorical(condition),
            },
            index=[f"spot_{i}" for i in range(n_obs)],
        ),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm["spatial"] = coords
    adata.layers["counts"] = counts.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    adata.var["highly_variable"] = True
    sc.pp.pca(adata, n_comps=10, random_state=42)
    sc.pp.neighbors(adata, n_neighbors=10, random_state=42)
    adata.write_h5ad(path)
    return path


class TestSpatialPlugin020(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.work = cls.tmpdir.name
        cls.h5ad = make_synthetic_spatial_h5ad(os.path.join(cls.work, "input.h5ad"))
        cls.gene_sets_path = os.path.join(cls.work, "gene_sets.json")
        with open(cls.gene_sets_path, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "domain0_prog": ["TGFB1", "VEGFA", "GENE00", "GENE01", "GENE02"],
                    "domain2_prog": ["CCR2", "KDR", "GENE16", "GENE17", "GENE18"],
                },
                handle,
            )

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def test_01_spatial_statistics_topology_and_ripley(self):
        out_h5 = os.path.join(self.work, "spatial_out.h5ad")
        figs = os.path.join(self.work, "figs_spatial")
        res = run_tool(
            "spatial_statistics.py",
            [
                "--input", self.h5ad,
                "--output", out_h5,
                "--figures_dir", figs,
                "--group_by", "leiden",
                "--n_perms", "50",
                "--cooc_interval", "10",
                "--ripley_mode", "L",
                "--ripley_simulations", "20",
                "--group_column", "condition",
            ],
        )
        self.assertIn("nhood_enrichment", res)
        self.assertIn("interaction_matrix", res)
        self.assertNotIn("skipped", res["interaction_matrix"])
        self.assertGreater(len(res["interaction_matrix"]["top_contact_pairs"]), 0)
        self.assertIn("centrality", res)
        self.assertNotIn("skipped", res["centrality"])
        self.assertIn("0", res["centrality"]["scores"])
        self.assertIn("ripley", res)
        self.assertNotIn("skipped", res["ripley"])
        self.assertEqual(res["ripley"]["mode"], "L")
        self.assertIn("0", res["ripley"]["per_cluster"])
        self.assertTrue(os.path.exists(out_h5))
        self.assertGreaterEqual(len(res["figures"]), 5)

    def test_02_niche_and_stratified_communication(self):
        niche_h5 = os.path.join(self.work, "niche_out.h5ad")
        figs_niche = os.path.join(self.work, "figs_niche")
        res_comp = run_tool(
            "niche.py",
            [
                "--input", self.h5ad,
                "--output", niche_h5,
                "--figures_dir", figs_niche,
                "--group_by", "leiden",
                "--niche_mode", "composition",
                "--niche_resolution", "0.5",
                "--group_column", "condition",
            ],
        )
        self.assertGreaterEqual(res_comp["niche"]["n_niches"], 2)
        self.assertIsNotNone(res_comp["niche"]["condition_dynamics"])
        self.assertEqual(res_comp["niche"]["condition_dynamics"]["group_column"], "condition")
        self.assertTrue(os.path.exists(niche_h5))

        # Expression mode smoke check.
        res_expr = run_tool(
            "niche.py",
            [
                "--input", self.h5ad,
                "--group_by", "leiden",
                "--niche_mode", "expression",
                "--niche_resolution", "0.5",
            ],
        )
        self.assertEqual(res_expr["niche"]["mode"], "expression")
        self.assertGreaterEqual(res_expr["niche"]["n_niches"], 2)

        # Niche-stratified cell-cell communication on niche_out.h5ad.
        ccc_h5 = os.path.join(self.work, "ccc_out.h5ad")
        figs_ccc = os.path.join(self.work, "figs_ccc")
        res_ccc = run_tool(
            "communication.py",
            [
                "--input", niche_h5,
                "--output", ccc_h5,
                "--figures_dir", figs_ccc,
                "--group_by", "leiden",
                "--stratify_by", "niche",
                "--n_perms", "50",
                "--ccc_threshold", "0.2",
            ],
        )
        self.assertEqual(res_ccc["communication"]["method"], "ligrec")
        self.assertEqual(res_ccc["communication"]["stratify_by"], "niche")
        self.assertIn("by_stratum", res_ccc["communication"])
        self.assertTrue(os.path.exists(ccc_h5))

    def test_03_trajectory_and_gene_scores(self):
        traj_h5 = os.path.join(self.work, "traj_out.h5ad")
        figs_traj = os.path.join(self.work, "figs_traj")
        res = run_tool(
            "trajectory.py",
            [
                "--input", self.h5ad,
                "--output", traj_h5,
                "--figures_dir", figs_traj,
                "--group_by", "leiden",
                "--root_cluster", "0",
                "--gene_sets_path", self.gene_sets_path,
                "--paga_threshold", "0.01",
            ],
        )
        self.assertEqual(res["trajectory"]["root_cluster"], "0")
        self.assertEqual(res["trajectory"]["root_selection"], "explicit")
        self.assertIn("0", res["trajectory"]["median_pseudotime_by_group"])
        self.assertEqual(res["gene_scores"]["programs_scored"], 2)
        self.assertTrue(os.path.exists(traj_h5))
        saved = ad.read_h5ad(traj_h5)
        self.assertIn("dpt_pseudotime", saved.obs.columns)
        self.assertIn("score_domain0_prog", saved.obs.columns)

    def test_04_notebook_render_and_self_check(self):
        nb_dir = os.path.join(self.work, "nb_run")
        os.makedirs(nb_dir, exist_ok=True)
        nb_path = os.path.join(nb_dir, "spatial_020.ipynb")
        res = run_tool(
            "build_notebook.py",
            [
                "--steps", "setup,load,spatial-neighbors,nhood,niche,trajectory,save",
                "--data", self.h5ad,
                "--output", nb_path,
                "--min_genes", "10",
                "--n_hvg", "40",
                "--n_pcs", "10",
                "--dims", "10",
                "--group_by", "leiden",
                "--n_perms", "30",
                "--cooc_interval", "10",
                "--ripley_simulations", "15",
                "--root_cluster", "0",
                "--self-check",
            ],
        )
        self.assertEqual(res["self_check"], "passed")
        self.assertTrue(os.path.exists(nb_path))
        self.assertIn("niche", res["steps"])
        self.assertIn("trajectory", res["steps"])


if __name__ == "__main__":
    unittest.main()
