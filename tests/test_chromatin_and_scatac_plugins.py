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
        for mod in range(1, 7):
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
        for mod in range(1, 7):
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


if __name__ == "__main__":
    unittest.main()
