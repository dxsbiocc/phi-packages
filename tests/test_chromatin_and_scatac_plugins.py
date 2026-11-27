#!/usr/bin/env python3
"""Regression tests for chromatin-profiling (0.1.0) and single-cell-atac (0.1.0):
- Fragment library integrity (steps.json <-> assets/fragments <-> references/step-notes)
- Dependency DAG validation and preset closure
- Full Python AST / R parser syntax verification across all scripts, fragments, and reports
- Live execution of scatac env_probe (+ --write_panel)
- Live execution of build_notebook.py across all 6 presets + --steps + --config for both plugins
  and syntax/structure verification of every rendered .ipynb
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHIP_SKILL = os.path.join(
    REPO_ROOT, "resources", "plugins", "chromatin-profiling", "skills", "chromatin-workflow"
)
SCATAC_SKILL = os.path.join(
    REPO_ROOT, "resources", "plugins", "single-cell-atac", "skills", "scatac-workflow"
)


def run_python(script_path, args):
    cmd = [sys.executable, script_path] + args
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if proc.returncode != 0:
        raise AssertionError(
            f"{os.path.basename(script_path)} failed (rc={proc.returncode})\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


def check_library_integrity(testcase, skill_dir, ext):
    steps_path = os.path.join(skill_dir, "assets", "steps.json")
    with open(steps_path, encoding="utf-8") as fh:
        reg = json.load(fh)

    step_ids = [s["id"] for s in reg["steps"]]
    testcase.assertEqual(len(step_ids), len(set(step_ids)), "duplicate step ids")
    step_set = set(step_ids)

    for step in reg["steps"]:
        sid = step["id"]
        frag = os.path.join(skill_dir, "assets", "fragments", f"{sid}.{ext}")
        note = os.path.join(skill_dir, "references", "step-notes", f"{sid}.md")
        testcase.assertTrue(os.path.exists(frag), f"missing fragment: {frag}")
        testcase.assertTrue(os.path.exists(note), f"missing step-note: {note}")
        for dep in step.get("requires", []):
            testcase.assertIn(dep, step_set, f"unknown dependency '{dep}' in step '{sid}'")

    for preset_id, preset in reg["presets"].items():
        for sid in preset["steps"]:
            testcase.assertIn(sid, step_set, f"preset {preset_id} references unknown step '{sid}'")

    return reg


class TestChromatinProfilingPlugin(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.work = cls.tmpdir.name
        # Minimal BED fixture for build_notebook.py --data
        cls.bed_path = os.path.join(cls.work, "sample_peaks.narrowPeak")
        with open(cls.bed_path, "w", encoding="utf-8") as fh:
            fh.write("chr1\t1000\t1500\tpeak1\t100\t.\t5.2\t10.1\t8.2\t250\n")

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def test_01_fragment_library_and_r_syntax(self):
        check_library_integrity(self, CHIP_SKILL, "R")

        # Verify every .py file compiles cleanly
        for root, _, files in os.walk(CHIP_SKILL):
            for fname in files:
                if fname.endswith(".py"):
                    path = os.path.join(root, fname)
                    with open(path, encoding="utf-8") as fh:
                        compile(fh.read(), path, "exec")

        # Verify every .R file parses cleanly if Rscript is available on host
        rscript = shutil.which("Rscript")
        if rscript:
            r_files = []
            for root, _, files in os.walk(CHIP_SKILL):
                for fname in sorted(files):
                    if fname.endswith(".R"):
                        r_files.append(os.path.join(root, fname))
            self.assertGreater(len(r_files), 15)
            expr = "; ".join(f'parse(file="{p}")' for p in r_files)
            proc = subprocess.run([rscript, "-e", expr], capture_output=True, text=True, timeout=60)
            self.assertEqual(proc.returncode, 0, f"R parse error:\n{proc.stderr}")

    def test_02_render_all_presets_and_config_override(self):
        builder = os.path.join(CHIP_SKILL, "scripts", "build_notebook.py")
        for mod in range(1, 8):
            out_nb = os.path.join(self.work, f"chip_mod{mod}.ipynb")
            res = run_python(builder, ["--module", str(mod), "--data", self.bed_path, "--output", out_nb])
            self.assertEqual(res["language"], "r")
            self.assertGreater(res["code_cells"], 2)
            with open(out_nb, encoding="utf-8") as fh:
                nb = json.load(fh)
            self.assertEqual(nb["nbformat"], 4)
            self.assertEqual(nb["metadata"]["kernelspec"]["name"], "plugin-chromatin-profiling-chip")

        # Custom --steps dependency closure + --config override
        cfg = os.path.join(self.work, "chip_cfg.json")
        with open(cfg, "w", encoding="utf-8") as fh:
            json.dump({"min_occupancy": 3, "fdr_threshold": 0.01}, fh)
        out_custom = os.path.join(self.work, "chip_custom.ipynb")
        res_c = run_python(
            builder,
            ["--steps", "consensus,save", "--data", self.bed_path, "--config", cfg, "--output", out_custom],
        )
        self.assertIn("setup", res_c["steps"])
        self.assertIn("consensus", res_c["steps"])
        with open(out_custom, encoding="utf-8") as fh:
            nb = json.load(fh)
        config_code = "".join(c for cell in nb["cells"] if cell.get("metadata", {}).get("phi_source") == "config" for c in cell["source"])
        self.assertIn("min_occupancy", config_code)
        self.assertIn("3", config_code)

    def test_03_live_chip_env_probe_and_consensus_when_diffbind_installed(self):
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript not installed")
        chk = subprocess.run(
            [
                rscript,
                "-e",
                'quit(status = if (requireNamespace("DiffBind", quietly=TRUE) && requireNamespace("ChIPseeker", quietly=TRUE)) 0 else 1)',
            ],
            capture_output=True,
            text=True,
        )
        if chk.returncode != 0:
            self.skipTest("DiffBind/ChIPseeker not installed in active R environment")

        probe = os.path.join(CHIP_SKILL, "scripts", "env_probe.R")
        proc = subprocess.run([rscript, probe], capture_output=True, text=True, timeout=120)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        res = json.loads(proc.stdout)
        self.assertGreaterEqual(res["packages_installed"], 6)
        self.assertTrue(any("seacr/callpeak" in g.get("next_step", "") for g in res.get("gaps", [])))

        # Run live chip_consensus on two synthetic euchromatic hg38 BED files (outside ENCODE telomere/centromere blacklist)
        pdir = os.path.join(self.work, "peaks_dir")
        os.makedirs(pdir, exist_ok=True)
        for sname in ("S1", "S2"):
            with open(os.path.join(pdir, f"{sname}.narrowPeak"), "w", encoding="utf-8") as fh:
                fh.write("chr1\t50000000\t50000500\tp1\t100\t.\t5.0\t10.0\t8.0\t250\n")
                fh.write("chr1\t50100000\t50100600\tp2\t100\t.\t6.0\t12.0\t9.0\t300\n")
        out_rds = os.path.join(self.work, "consensus.rds")
        out_bed = os.path.join(self.work, "consensus.bed")
        figs_dir = os.path.join(self.work, "figs_chip")
        cons_script = os.path.join(CHIP_SKILL, "scripts", "consensus_peaks.R")
        cproc = subprocess.run(
            [
                rscript,
                cons_script,
                "--input",
                pdir,
                "--genome",
                "hg38",
                "--min_occupancy",
                "2",
                "--output",
                out_rds,
                "--figures_dir",
                figs_dir,
            ],
            capture_output=True,
            text=True,
            timeout=180,
        )
        self.assertEqual(cproc.returncode, 0, f"STDOUT:\n{cproc.stdout}\nSTDERR:\n{cproc.stderr}")
        cres = json.loads(cproc.stdout)
        self.assertEqual(cres["n_consensus"], 2)
        self.assertTrue(os.path.isfile(out_rds))
        self.assertTrue(os.path.isfile(out_bed))
        self.assertGreaterEqual(len(cres.get("artifacts", [])), 1)
        for art in cres.get("artifacts", []):
            self.assertTrue(os.path.isfile(art + ".phi-artifact.json"))

        # Run live chip_super_enhancer on bimodal stitched enhancers + typical enhancers
        se_bed = os.path.join(self.work, "enhancers.narrowPeak")
        with open(se_bed, "w", encoding="utf-8") as fh:
            for i in range(12):
                pos = 1000000 + i * 50000
                fh.write(f"chr2\t{pos}\t{pos + 500}\tte_{i}\t10\t.\t{2.0 + 0.1 * i}\t5.0\t3.0\t250\n")
            for j in range(6):
                pos = 5000000 + j * 2000
                fh.write(f"chr2\t{pos}\t{pos + 600}\tse_const_{j}\t200\t.\t85.0\t20.0\t15.0\t300\n")
        se_tsv = os.path.join(self.work, "super_enhancers.tsv")
        se_script = os.path.join(CHIP_SKILL, "scripts", "super_enhancer.R")
        se_proc = subprocess.run(
            [
                rscript,
                se_script,
                "--input",
                se_bed,
                "--se_output",
                se_tsv,
                "--figures_dir",
                os.path.join(self.work, "figs_se"),
                "--stitch_distance",
                "12500",
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertEqual(se_proc.returncode, 0, f"STDOUT:\n{se_proc.stdout}\nSTDERR:\n{se_proc.stderr}")
        se_res = json.loads(se_proc.stdout)
        self.assertEqual(se_res["n_super_enhancers"], 1)
        self.assertGreaterEqual(se_res["n_typical_enhancers"], 10)
        self.assertTrue(os.path.isfile(se_tsv))


class TestSingleCellAtacPlugin(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.work = cls.tmpdir.name
        cls.frag_path = os.path.join(cls.work, "fragments.tsv.gz")
        with open(cls.frag_path, "wb") as fh:
            fh.write(b"")  # build_notebook only checks existence & fingerprint without --self-check

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def test_01_fragment_library_and_python_syntax(self):
        check_library_integrity(self, SCATAC_SKILL, "py")

        py_files = []
        for root, _, files in os.walk(SCATAC_SKILL):
            for fname in sorted(files):
                if fname.endswith(".py"):
                    path = os.path.join(root, fname)
                    py_files.append(path)
                    with open(path, encoding="utf-8") as fh:
                        compile(fh.read(), path, "exec")
        self.assertGreater(len(py_files), 20)

    def test_02_env_probe_and_starter_panel(self):
        probe = os.path.join(SCATAC_SKILL, "scripts", "env_probe.py")
        panel_out = os.path.join(self.work, "starter_panel.json")
        res = run_python(probe, ["--write_panel", panel_out])
        self.assertIn("packages", res)
        self.assertIn("modules", res)
        self.assertEqual(len(res["modules"]), 6)
        self.assertTrue(os.path.exists(panel_out))
        with open(panel_out, encoding="utf-8") as fh:
            panel = json.load(fh)
        self.assertIsInstance(panel, dict)
        self.assertGreaterEqual(len(panel), 3)

    def test_03_render_all_presets_and_compile_every_cell(self):
        builder = os.path.join(SCATAC_SKILL, "scripts", "build_notebook.py")
        for mod in range(1, 8):
            out_nb = os.path.join(self.work, f"scatac_mod{mod}.ipynb")
            res = run_python(builder, ["--module", str(mod), "--data", self.frag_path, "--output", out_nb])
            self.assertEqual(res["language"], "python")
            self.assertGreater(res["code_cells"], 2)
            with open(out_nb, encoding="utf-8") as fh:
                nb = json.load(fh)
            self.assertEqual(nb["nbformat"], 4)
            self.assertEqual(nb["metadata"]["kernelspec"]["name"], "plugin-single-cell-atac-scatac")
            # Compile every rendered Python code cell to catch broken template expansion
            for idx, cell in enumerate(nb["cells"]):
                if cell["cell_type"] == "code":
                    src = "".join(cell["source"])
                    compile(src, f"scatac_mod{mod}_cell_{idx}", "exec")

    def test_04_scatac_tf_deviation_and_peak2gene_end_to_end(self):
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(42)
        n_peaks, n_cells = 40, 24
        peaks = [f"chr1:{100000 + i * 10000}-{100500 + i * 10000}" for i in range(30)] + [
            f"chr2:{200000 + j * 10000}-{200500 + j * 10000}" for j in range(10)
        ]
        cells = [f"C{k:02d}" for k in range(n_cells)]
        # Plant active TF1 peaks 0..7 in cells 0..11, and reciprocal ClusterB peaks 8..15 in cells 12..23
        X = rng.poisson(2.0, size=(n_peaks, n_cells)).astype(float)
        X[:8, :12] += 18.0
        X[8:16, 12:] += 18.0
        peak_df = pd.DataFrame(X, index=peaks, columns=cells)
        peak_tsv = os.path.join(self.work, "peaks_counts.tsv")
        peak_df.to_csv(peak_tsv, sep="\t")

        meta_df = pd.DataFrame(
            {"leiden": ["ClusterA" if k < 12 else "ClusterB" for k in range(n_cells)]},
            index=cells,
        )
        meta_tsv = os.path.join(self.work, "atac_meta.tsv")
        meta_df.to_csv(meta_tsv, sep="\t")

        motif_df = pd.DataFrame(
            {
                "TF_ACTIVE": [1 if i < 8 else 0 for i in range(n_peaks)],
                "TF_NOISE": [1 if 20 <= i < 28 else 0 for i in range(n_peaks)],
            },
            index=peaks,
        )
        motif_tsv = os.path.join(self.work, "motif_matches.tsv")
        motif_df.to_csv(motif_tsv, sep="\t")

        gc_df = pd.DataFrame({"peak": peaks, "gc": np.linspace(0.35, 0.65, n_peaks)})
        gc_tsv = os.path.join(self.work, "peak_gc.tsv")
        gc_df.to_csv(gc_tsv, sep="\t", index=False)

        tfdev_tsv = os.path.join(self.work, "tf_dev.tsv")
        tf_res = run_python(
            os.path.join(SCATAC_SKILL, "scripts", "tf_deviation.py"),
            [
                "--input", peak_tsv,
                "--meta_path", meta_tsv,
                "--group_by", "leiden",
                "--motif_matches", motif_tsv,
                "--peak_gc", gc_tsv,
                "--tfdev_output", tfdev_tsv,
                "--figures_dir", os.path.join(self.work, "figs_tfdev"),
                "--n_background", "30",
            ],
        )
        self.assertEqual(tf_res["method"], "chromvar_gc_matched")
        self.assertEqual(tf_res["gc_matching"], "peak_gc_table")
        dev_mat = pd.read_csv(tfdev_tsv, sep="\t", index_col=0)
        self.assertGreater(float(dev_mat.loc[cells[:12], "TF_ACTIVE"].mean()), 2.0)
        self.assertTrue(os.path.isfile(tfdev_tsv + ".phi-artifact.json"))

        # Plant cis enhancer link between peak 0 (chr1:100000-100500) and GENE_ cis on chr1 (TSS=150000)
        rna_mat = rng.normal(5.0, 0.4, size=(4, n_cells))
        rna_mat[0, :] = X[0, :] * 0.8 + rng.normal(0.0, 0.3, size=n_cells)
        genes = ["GENE_TARGET", "GENE_NULL1", "GENE_NULL2", "GENE_CHR2"]
        rna_tsv = os.path.join(self.work, "rna_counts.tsv")
        pd.DataFrame(rna_mat, index=genes, columns=cells).to_csv(rna_tsv, sep="\t")

        tss_tsv = os.path.join(self.work, "gene_tss.tsv")
        pd.DataFrame({
            "gene": genes,
            "chrom": ["chr1", "chr1", "chr1", "chr2"],
            "tss": [150000, 120000, 180000, 220000],
        }).to_csv(tss_tsv, sep="\t", index=False)

        p2g_tsv = os.path.join(self.work, "p2g_links.tsv")
        p2g_res = run_python(
            os.path.join(SCATAC_SKILL, "scripts", "peak2gene.py"),
            [
                "--input", peak_tsv,
                "--rna_input", rna_tsv,
                "--gene_tss", tss_tsv,
                "--max_distance", "100000",
                "--n_null", "200",
                "--links_output", p2g_tsv,
                "--figures_dir", os.path.join(self.work, "figs_p2g"),
            ],
        )
        self.assertEqual(p2g_res["pairing_mode"], "paired_barcodes")
        links_df = pd.read_csv(p2g_tsv, sep="\t")
        top_hit = links_df[(links_df["gene"] == "GENE_TARGET") & (links_df["peak"] == peaks[0])].iloc[0]
        self.assertGreater(float(top_hit["pearson_r"]), 0.7)
        self.assertLess(float(top_hit["fdr"]), 0.05)


if __name__ == "__main__":
    unittest.main()
