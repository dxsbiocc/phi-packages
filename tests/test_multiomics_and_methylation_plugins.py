#!/usr/bin/env python3
"""Contract, layout, AST, and end-to-end tests for multi-omics (0.2.0) and methylation (0.1.1)."""

from __future__ import annotations

import gzip
import importlib.util
import json
import py_compile
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MO_DIR = ROOT / "resources/plugins/multi-omics"
MO_SCRIPTS = MO_DIR / "skills/multiomics-workflow/scripts"
MO_SCHEMAS = MO_DIR / "skills/multiomics-workflow/schemas"

METH_DIR = ROOT / "resources/plugins/methylation"
METH_SCRIPTS = METH_DIR / "skills/methylation-workflow/scripts"
METH_SCHEMAS = METH_DIR / "skills/methylation-workflow/schemas"

HAS_SCIPY_STACK = all(
    importlib.util.find_spec(pkg) is not None for pkg in ("numpy", "pandas", "scipy", "matplotlib")
)


def run_json(script: Path, args: list[str], cwd: Path) -> dict:
    proc = subprocess.run(
        [sys.executable, "-B", str(script), *args],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"{script.name} failed (rc={proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


class TestMultiOmicsAndMethylationContracts(unittest.TestCase):
    def test_multiomics_layout_locks_schemas_and_ast(self):
        manifest = (MO_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("id: multi-omics", manifest)
        self.assertIn("version: 0.2.0", manifest)
        self.assertIn("toolPrefix: mo", manifest)
        for req in ("icon.svg", "README.md", "CHANGELOG.md", "agents/MultiOmics.md"):
            self.assertTrue((MO_DIR / req).is_file(), f"missing {req}")
        for plat in ("darwin-arm64", "darwin-x64", "linux-x64"):
            lock = MO_DIR / f"environments/multiomics/locks/{plat}.txt"
            self.assertTrue(lock.is_file())
            self.assertIn("@EXPLICIT", lock.read_text(encoding="utf-8"))

        self.assertEqual(
            {p.name for p in MO_SCHEMAS.glob("*.json")},
            {
                "inspect-result.json",
                "rna-protein-result.json",
                "pair-assoc-result.json",
                "factor-integration-result.json",
                "subtype-result.json",
            },
        )
        for sp in MO_SCHEMAS.glob("*.json"):
            self.assertEqual(json.loads(sp.read_text(encoding="utf-8"))["type"], "object")
        for py in MO_SCRIPTS.glob("*.py"):
            py_compile.compile(str(py), doraise=True)

    def test_methylation_layout_locks_schemas_and_ast(self):
        manifest = (METH_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("id: methylation", manifest)
        self.assertIn("version: 0.1.1", manifest)
        self.assertIn("toolPrefix: meth", manifest)
        for req in ("icon.svg", "README.md", "CHANGELOG.md", "agents/Methylation.md"):
            self.assertTrue((METH_DIR / req).is_file(), f"missing {req}")
        for plat in ("darwin-arm64", "darwin-x64", "linux-x64"):
            lock = METH_DIR / f"environments/methylation/locks/{plat}.txt"
            self.assertTrue(lock.is_file())
            self.assertIn("@EXPLICIT", lock.read_text(encoding="utf-8"))

        self.assertEqual(
            {p.name for p in METH_SCHEMAS.glob("*.json")},
            {"inspect-result.json", "merge-result.json", "qc-result.json", "dmr-result.json", "annotate-result.json"},
        )
        for sp in METH_SCHEMAS.glob("*.json"):
            self.assertEqual(json.loads(sp.read_text(encoding="utf-8"))["type"], "object")
        for py in METH_SCRIPTS.glob("*.py"):
            py_compile.compile(str(py), doraise=True)


@unittest.skipUnless(HAS_SCIPY_STACK, "numpy/pandas/scipy/matplotlib required for live E2E")
class TestMultiOmicsAndMethylationE2E(unittest.TestCase):
    def test_multiomics_full_chain(self):
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(101)
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            n = 24
            samples = [f"S{i:02d}" for i in range(1, n + 1)]
            grp = ["Case"] * 12 + ["Control"] * 12
            age = np.round(rng.normal(55, 8, size=n), 1)
            meta = pd.DataFrame({"sample": samples, "group": grp, "age": age})
            meta_p = tdir / "meta.tsv"
            meta.to_csv(meta_p, sep="\t", index=False)

            # Shared latent factor driving RNA, Protein, and Metabolome
            latent = np.concatenate([np.ones(12), -np.ones(12)]) + rng.normal(0, 0.35, size=n)
            genes = [f"GENE{i:02d}" for i in range(1, 21)]
            rna_mat = pd.DataFrame(
                np.outer(np.linspace(1.2, -1.0, 20), latent) + rng.normal(0, 0.4, size=(20, n)) + 8.0,
                index=genes,
                columns=samples,
            )
            prot_mat = pd.DataFrame(
                np.outer(np.linspace(1.0, -0.8, 20), latent) + rng.normal(0, 0.4, size=(20, n)) + 20.0,
                index=genes,
                columns=samples,
            )
            # Compositional microbiome (sums to 100)
            raw_taxa = np.exp(np.outer(np.linspace(0.9, -0.9, 12), latent) + rng.normal(0, 0.3, size=(12, n)))
            taxa_mat = pd.DataFrame(
                100.0 * raw_taxa / raw_taxa.sum(axis=0, keepdims=True),
                index=[f"Taxon_{i:02d}" for i in range(1, 13)],
                columns=samples,
            )
            rna_p = tdir / "rna.tsv"
            prot_p = tdir / "prot.tsv"
            taxa_p = tdir / "taxa.tsv"
            rna_mat.to_csv(rna_p, sep="\t")
            prot_mat.to_csv(prot_p, sep="\t")
            taxa_mat.to_csv(taxa_p, sep="\t")

            # 1) mo_inspect
            ins = run_json(
                MO_SCRIPTS / "mo_inspect.py",
                [
                    "--matrices", f"rna={rna_p}",
                    "--matrices", f"protein={prot_p}",
                    "--matrices", f"microbiome={taxa_p}",
                    "--metadata", str(meta_p),
                    "--group_col", "group",
                ],
                cwd=tdir,
            )
            self.assertEqual(ins["n_shared_samples"], n)
            scales = {m["modality"]: m["inferred_scale"] for m in ins["modalities"]}
            self.assertEqual(scales["microbiome"], "relative_compositional")

            # 2) mo_rna_protein (DE 9-quadrant + sample-level correlation)
            rna_de = pd.DataFrame({
                "gene": genes,
                "log2FoldChange": np.concatenate([[2.2, 1.8, 2.0, -2.1], rng.normal(0, 0.3, size=16)]),
                "padj": np.concatenate([[1e-5, 1e-4, 1e-4, 2e-5], np.full(16, 0.5)]),
            })
            prot_de = pd.DataFrame({
                "gene": genes,
                "log2FC": np.concatenate([[1.5, 0.05, -1.2, -1.4], rng.normal(0, 0.2, size=16)]),
                "fdr": np.concatenate([[1e-4, 0.7, 1e-3, 1e-4], np.full(16, 0.6)]),
            })
            rde_p = tdir / "rna_de.tsv"
            pde_p = tdir / "prot_de.tsv"
            rna_de.to_csv(rde_p, sep="\t", index=False)
            prot_de.to_csv(pde_p, sep="\t", index=False)

            rp = run_json(
                MO_SCRIPTS / "mo_rna_protein.py",
                [
                    "--rna_de", str(rde_p),
                    "--prot_de", str(pde_p),
                    "--rna_matrix", str(rna_p),
                    "--prot_matrix", str(prot_p),
                    "--output_dir", str(tdir / "rp_out"),
                ],
                cwd=tdir,
            )
            q_cnt = rp["de_concordance"]["quadrant_counts"]
            self.assertIn("concordant_up", q_cnt)
            self.assertIn("rna_only_up_buffered", q_cnt)
            self.assertIsNotNone(rp["sample_level_correlation"])
            for art in rp["artifacts"]:
                self.assertTrue(Path(art + ".phi-artifact.json").is_file())

            # 3) mo_pair_assoc (CLR microbiome vs protein with age covariate)
            pa = run_json(
                MO_SCRIPTS / "mo_pair_assoc.py",
                [
                    "--matrix_x", str(taxa_p),
                    "--matrix_y", str(prot_p),
                    "--name_x", "Microbiome",
                    "--name_y", "Proteome",
                    "--transform_x", "clr",
                    "--metadata", str(meta_p),
                    "--group_col", "group",
                    "--covariates", "age",
                    "--permutations", "99",
                    "--output_dir", str(tdir / "pa_out"),
                ],
                cwd=tdir,
            )
            self.assertEqual(pa["network"]["method"], "partial_spearman")
            self.assertLess(pa["procrustes"]["pvalue"], 0.05)

            # 4) mo_factor_integration (3 modalities)
            fi = run_json(
                MO_SCRIPTS / "mo_factor_integration.py",
                [
                    "--matrices", f"rna={rna_p}",
                    "--matrices", f"protein={prot_p}",
                    "--matrices", f"microbiome={taxa_p}",
                    "--metadata", str(meta_p),
                    "--group_col", "group",
                    "--n_factors", "3",
                    "--output_dir", str(tdir / "fi_out"),
                ],
                cwd=tdir,
            )
            self.assertEqual(fi["n_factors"], 3)
            self.assertEqual(fi["factor_variance_explained"][0]["kind"], "shared_multiomics")

            # 5) mo_subtype (SNF multi-omics molecular subtyping)
            snf = run_json(
                MO_SCRIPTS / "mo_subtype.py",
                [
                    "--matrices", f"rna={rna_p}",
                    "--matrices", f"protein={prot_p}",
                    "--matrices", f"microbiome={taxa_p}",
                    "--metadata", str(meta_p),
                    "--group_col", "group",
                    "--n_neighbors", "8",
                    "--n_iters", "12",
                    "--output_dir", str(tdir / "snf_out"),
                ],
                cwd=tdir,
            )
            self.assertEqual(snf["tool"], "mo_subtype")
            self.assertEqual(snf["n_shared_samples"], n)
            self.assertEqual(snf["n_modalities"], 3)
            self.assertGreater(snf["silhouette_score"], 0.2)
            self.assertGreater(snf["group_association"]["ari"], 0.75)
            for art in snf["artifacts"]:
                self.assertTrue(Path(art + ".phi-artifact.json").is_file())

    def test_methylation_full_chain(self):
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            cov_dir = tdir / "bismark"
            cov_dir.mkdir()
            samples = ["Tumor_1", "Tumor_2", "Tumor_3", "Norm_1", "Norm_2", "Norm_3"]
            # 4 adjacent CpGs in TP53 promoter CpG island (hypermethylated in Tumor) + 16 background CpGs
            dmr_pos = [1050, 1120, 1190, 1260]
            bg_pos = [5000 + i * 800 for i in range(16)]

            for s_idx, sname in enumerate(samples):
                is_tumor = sname.startswith("Tumor")
                lines = []
                for p in dmr_pos:
                    m_cnt = (42 + s_idx) if is_tumor else (5 + (s_idx % 2))
                    u_cnt = 50 - m_cnt
                    pct = 100.0 * m_cnt / 50.0
                    lines.append(f"chr17\t{p}\t{p}\t{pct:.1f}\t{m_cnt}\t{u_cnt}")
                for b_i, p in enumerate(bg_pos):
                    m_cnt = 10 + ((b_i + s_idx) % 5)
                    u_cnt = 40 - m_cnt
                    pct = 100.0 * m_cnt / 40.0
                    lines.append(f"chr17\t{p}\t{p}\t{pct:.1f}\t{m_cnt}\t{u_cnt}")
                with gzip.open(cov_dir / f"{sname}.bismark.cov.gz", "wt", encoding="utf-8") as fh:
                    fh.write("\n".join(lines) + "\n")

            meta = pd.DataFrame({
                "sample": samples,
                "condition": ["Tumor", "Tumor", "Tumor", "Normal", "Normal", "Normal"],
                "age": [61, 65, 59, 60, 63, 58],
            })
            meta_p = tdir / "meta.tsv"
            meta.to_csv(meta_p, sep="\t", index=False)

            # 1) meth_inspect on directory
            ins = run_json(
                METH_SCRIPTS / "meth_inspect.py",
                ["--input", str(cov_dir), "--metadata", str(meta_p), "--group_col", "condition"],
                cwd=tdir,
            )
            self.assertEqual(ins["n_samples"], 6)

            # 2) meth_merge
            mrg_dir = tdir / "merged"
            mrg = run_json(
                METH_SCRIPTS / "meth_merge.py",
                ["--inputs", str(cov_dir), "--min_cov", "10", "--output_dir", str(mrg_dir)],
                cwd=tdir,
            )
            self.assertEqual(mrg["n_cpgs_retained"], 20)
            beta_p = mrg["outputs"]["beta_matrix"]

            # 3) meth_qc
            qc = run_json(
                METH_SCRIPTS / "meth_qc.py",
                [
                    "--matrix", beta_p,
                    "--coverage_matrix", mrg["outputs"]["coverage_matrix"],
                    "--metadata", str(meta_p),
                    "--group_col", "condition",
                    "--output_dir", str(tdir / "qc"),
                ],
                cwd=tdir,
            )
            self.assertEqual(qc["input_scale"], "beta")

            # 4) meth_dmr (with age covariate)
            dmr = run_json(
                METH_SCRIPTS / "meth_dmr.py",
                [
                    "--matrix", beta_p,
                    "--metadata", str(meta_p),
                    "--group_col", "condition",
                    "--case", "Tumor",
                    "--control", "Normal",
                    "--covariates", "age",
                    "--delta_beta", "0.20",
                    "--min_cpgs", "3",
                    "--output_dir", str(tdir / "dmr"),
                ],
                cwd=tdir,
            )
            self.assertEqual(dmr["dmc_summary"]["n_hyper_dmcs"], 4)
            self.assertEqual(dmr["dmr_summary"]["n_hyper_dmrs"], 1)

            # 5) meth_annotate (GTF + CpG Island BED)
            gtf = tdir / "ref.gtf"
            gtf.write_text(
                'chr17\tENSEMBL\tgene\t1200\t8000\t.\t+\t.\tgene_id "ENSG1"; gene_name "TP53";\n'
                'chr17\tENSEMBL\tgene\t9000\t15000\t.\t+\t.\tgene_id "ENSG2"; gene_name "WRAP53";\n',
                encoding="utf-8",
            )
            cgi = tdir / "cgi.bed"
            cgi.write_text("chr17\t900\t1500\tCGI_1\n", encoding="utf-8")

            ann = run_json(
                METH_SCRIPTS / "meth_annotate.py",
                [
                    "--input", dmr["outputs"]["dmr_table"],
                    "--gtf", str(gtf),
                    "--cpg_islands", str(cgi),
                    "--output_dir", str(tdir / "annot"),
                ],
                cwd=tdir,
            )
            self.assertEqual(ann["genomic_feature_counts"].get("Promoter"), 1)
            self.assertEqual(ann["cpg_island_context_counts"].get("Island"), 1)


if __name__ == "__main__":
    unittest.main()
