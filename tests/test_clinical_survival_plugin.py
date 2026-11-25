#!/usr/bin/env python3
"""Contract, layout, AST, and end-to-end tests for the clinical-survival plugin (0.1.0)."""

from __future__ import annotations

import importlib.util
import json
import py_compile
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLUGIN_DIR = ROOT / "resources/plugins/clinical-survival"
SKILL_DIR = PLUGIN_DIR / "skills/survival-workflow"
SCRIPTS_DIR = SKILL_DIR / "scripts"
SCHEMAS_DIR = SKILL_DIR / "schemas"

HAS_SCIPY_STACK = all(
    importlib.util.find_spec(pkg) is not None for pkg in ("numpy", "pandas", "scipy", "matplotlib")
)


def run_json(script_name: str, args: list[str], cwd: Path) -> dict:
    script = SCRIPTS_DIR / script_name
    proc = subprocess.run(
        [sys.executable, str(script), *args],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"{script_name} failed (rc={proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


class TestClinicalSurvivalPluginContract(unittest.TestCase):
    def test_plugin_layout_locks_schemas_and_ast(self):
        manifest = (PLUGIN_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("id: clinical-survival", manifest)
        self.assertIn("version: 0.1.0", manifest)
        self.assertIn("toolPrefix: clin", manifest)

        for req in ("icon.svg", "README.md", "CHANGELOG.md", "agents/ClinicalSurvival.md"):
            self.assertTrue((PLUGIN_DIR / req).is_file(), f"missing {req}")

        for plat in ("darwin-arm64", "darwin-x64", "linux-x64"):
            lock = PLUGIN_DIR / f"environments/clinical/locks/{plat}.txt"
            self.assertTrue(lock.is_file(), f"missing lock {plat}")
            self.assertIn("@EXPLICIT", lock.read_text(encoding="utf-8"))

        expected_schemas = {
            "inspect-result.json",
            "table1-result.json",
            "km-result.json",
            "cox-result.json",
            "time-roc-result.json",
        }
        actual_schemas = {p.name for p in SCHEMAS_DIR.glob("*.json")}
        self.assertEqual(actual_schemas, expected_schemas)
        for sp in SCHEMAS_DIR.glob("*.json"):
            data = json.loads(sp.read_text(encoding="utf-8"))
            self.assertEqual(data["type"], "object")

        expected_scripts = {
            "_common.py",
            "clin_inspect.py",
            "clin_table1.py",
            "clin_km.py",
            "clin_cox.py",
            "clin_time_roc.py",
        }
        actual_scripts = {p.name for p in SCRIPTS_DIR.glob("*.py")}
        self.assertEqual(actual_scripts, expected_scripts)
        for py in SCRIPTS_DIR.glob("*.py"):
            py_compile.compile(str(py), doraise=True)


@unittest.skipUnless(HAS_SCIPY_STACK, "numpy/pandas/scipy/matplotlib required for live survival E2E")
class TestClinicalSurvivalWorkflowE2E(unittest.TestCase):
    def test_full_clinical_survival_chain(self):
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(2026)
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            n = 90
            samples = [f"PT_{i:03d}" for i in range(1, n + 1)]
            arm = np.array(["Experimental"] * 45 + ["Control"] * 45)
            age = np.round(rng.normal(62.0, 9.5, size=n), 1)
            stage = rng.choice(["Stage_I_II", "Stage_III_IV"], size=n, p=[0.55, 0.45])
            gene_a = np.round(rng.normal(0.0, 1.0, size=n) + np.where(arm == "Control", 0.8, -0.4), 3)
            gene_b = np.round(rng.normal(0.0, 1.0, size=n), 3)

            # True log-hazard driven by arm, stage, and gene_a
            linpred = (
                0.75 * (arm == "Control").astype(float)
                + 0.65 * (stage == "Stage_III_IV").astype(float)
                + 0.80 * gene_a
                + 0.02 * (age - 60.0)
            )
            true_t = rng.exponential(scale=36.0 * np.exp(-linpred))
            cens_t = rng.uniform(12.0, 72.0, size=n)
            obs_t = np.round(np.maximum(0.5, np.minimum(true_t, cens_t)), 2)
            event_num = (true_t <= cens_t).astype(int)
            os_status = np.where(event_num == 1, "1:DECEASED", "0:LIVING")

            clin_df = pd.DataFrame(
                {
                    "patient_id": samples,
                    "os_months": obs_t,
                    "os_status": os_status,
                    "arm": arm,
                    "age": age,
                    "stage": stage,
                    "precomputed_risk": np.round(linpred, 4),
                }
            )
            clin_path = tdir / "clinical.tsv"
            clin_df.to_csv(clin_path, sep="\t", index=False)

            mat_df = pd.DataFrame(
                [gene_a, gene_b],
                index=["GENE_A", "GENE_B"],
                columns=samples,
            )
            mat_df.index.name = "gene_symbol"
            mat_path = tdir / "expr.tsv"
            mat_df.to_csv(mat_path, sep="\t")

            # 1) clin_inspect
            ins = run_json(
                "clin_inspect.py",
                [
                    "--clinical",
                    str(clin_path),
                    "--time_col",
                    "os_months",
                    "--event_col",
                    "os_status",
                    "--matrix",
                    str(mat_path),
                ],
                cwd=tdir,
            )
            self.assertEqual(ins["n_valid_survival"], n)
            self.assertGreater(ins["n_events"], 25)
            self.assertEqual(ins["inferred_time_unit"], "months")
            self.assertEqual(ins["matrix_overlap"]["n_matched_clinical_samples"], n)

            # 2) clin_table1 (using Phi engine style repeated --variables flags)
            t1_dir = tdir / "table1"
            t1 = run_json(
                "clin_table1.py",
                [
                    "--clinical",
                    str(clin_path),
                    "--group_col",
                    "arm",
                    "--variables",
                    "age",
                    "--variables",
                    "stage",
                    "--output_dir",
                    str(t1_dir),
                ],
                cwd=tdir,
            )
            self.assertEqual(len(t1["variables"]), 2)
            for art in t1["artifacts"]:
                self.assertTrue(Path(art + ".phi-artifact.json").is_file())

            # 3) clin_km (biomarker median split from molecular matrix + repeated --landmarks)
            km_dir = tdir / "km"
            km = run_json(
                "clin_km.py",
                [
                    "--clinical",
                    str(clin_path),
                    "--time_col",
                    "os_months",
                    "--event_col",
                    "os_status",
                    "--biomarker_col",
                    "GENE_A",
                    "--matrix",
                    str(mat_path),
                    "--landmarks",
                    "12",
                    "--landmarks",
                    "36",
                    "--output_dir",
                    str(km_dir),
                ],
                cwd=tdir,
            )
            self.assertLess(km["logrank_pvalue"], 0.05)
            self.assertEqual(len(km["landmark_survival"]), 4)
            for art in km["artifacts"]:
                self.assertTrue(Path(art + ".phi-artifact.json").is_file())

            # 4) clin_cox (both univariable + multivariable + Schoenfeld PH + boolean flag without value)
            cox_dir = tdir / "cox"
            cox = run_json(
                "clin_cox.py",
                [
                    "--clinical",
                    str(clin_path),
                    "--time_col",
                    "os_months",
                    "--event_col",
                    "os_status",
                    "--covariates",
                    "arm",
                    "--covariates",
                    "stage",
                    "--covariates",
                    "GENE_A",
                    "--matrix",
                    str(mat_path),
                    "--standardize_continuous",
                    "--output_dir",
                    str(cox_dir),
                ],
                cwd=tdir,
            )
            self.assertTrue(cox["multivariable"]["converged"])
            self.assertGreater(cox["multivariable"]["c_index"], 0.65)
            self.assertIsNotNone(cox["multivariable"]["global_ph_pvalue"])
            for art in cox["artifacts"]:
                self.assertTrue(Path(art + ".phi-artifact.json").is_file())

            # 5) clin_time_roc (leak-free K-fold CV L2-Cox signature)
            roc_dir = tdir / "roc"
            roc = run_json(
                "clin_time_roc.py",
                [
                    "--clinical",
                    str(clin_path),
                    "--time_col",
                    "os_months",
                    "--event_col",
                    "os_status",
                    "--feature_cols",
                    "age",
                    "--feature_cols",
                    "GENE_A",
                    "--feature_cols",
                    "GENE_B",
                    "--matrix",
                    str(mat_path),
                    "--landmarks",
                    "12",
                    "--landmarks",
                    "24",
                    "--landmarks",
                    "36",
                    "--output_dir",
                    str(roc_dir),
                ],
                cwd=tdir,
            )
            self.assertEqual(roc["evaluation_mode"], "oof_cv_l2_cox")
            self.assertGreater(roc["c_index"], 0.60)
            self.assertEqual(len(roc["landmarks"]), 3)
            for art in roc["artifacts"]:
                self.assertTrue(Path(art + ".phi-artifact.json").is_file())


if __name__ == "__main__":
    unittest.main()
