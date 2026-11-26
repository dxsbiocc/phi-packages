#!/usr/bin/env python3
"""Contract, layout, AST, schema, and end-to-end tests for the unified clinical plugin (0.1.0)."""

from __future__ import annotations

import ast
import importlib.util
import json
import math
import py_compile
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLUGIN_DIR = ROOT / "resources/plugins/clinical"
SURV_DIR = PLUGIN_DIR / "skills/survival-workflow"
PRED_DIR = PLUGIN_DIR / "skills/prediction-causal-workflow"

HAS_SCIPY_STACK = all(
    importlib.util.find_spec(pkg) is not None
    for pkg in ("numpy", "pandas", "scipy", "sklearn", "matplotlib")
)
HAS_JSONSCHEMA = importlib.util.find_spec("jsonschema") is not None


def run_tool(skill_dir: Path, script_name: str, args: list[str], cwd: Path) -> dict:
    script = skill_dir / "scripts" / script_name
    proc = subprocess.run(
        [sys.executable, "-B", str(script), *args],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"{script_name} failed (rc={proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise AssertionError(f"{script_name} did not emit single valid JSON: {exc}\nSTDOUT:\n{proc.stdout}") from exc


def validate_output_schema(payload: dict, schema_path: Path) -> None:
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    if HAS_JSONSCHEMA:
        import jsonschema
        jsonschema.validate(instance=payload, schema=schema)
    for req in schema.get("required", []):
        if req not in payload:
            raise AssertionError(f"missing required key '{req}' for schema {schema_path.name}")
    props = schema.get("properties", {})
    type_map = {
        "string": str,
        "integer": int,
        "number": (int, float),
        "boolean": bool,
        "array": list,
        "object": dict,
        "null": type(None),
    }
    for key, val in payload.items():
        if key not in props:
            continue
        spec = props[key]
        t_decl = spec.get("type")
        if t_decl:
            allowed = t_decl if isinstance(t_decl, list) else [t_decl]
            allowed_types = tuple(type_map[t] for t in allowed)
            if not isinstance(val, allowed_types) or (isinstance(val, bool) and "boolean" not in allowed):
                raise AssertionError(f"{schema_path.name}: key '{key}' got {type(val).__name__}, expected {allowed}")
        if "enum" in spec and val not in spec["enum"]:
            raise AssertionError(f"{schema_path.name}: key '{key}'={val!r} not in enum {spec['enum']}")


def assert_artifacts_valid(testcase: unittest.TestCase, payload: dict, expected_skill: str) -> None:
    testcase.assertGreater(len(payload["artifacts"]), 0)
    for art in payload["artifacts"]:
        art_p = Path(art)
        sidecar_p = Path(art + ".phi-artifact.json")
        testcase.assertTrue(art_p.is_file(), f"missing artifact file: {art}")
        testcase.assertTrue(sidecar_p.is_file(), f"missing artifact sidecar: {sidecar_p}")
        desc = json.loads(sidecar_p.read_text(encoding="utf-8"))
        testcase.assertEqual(desc["contractVersion"], "1.0.0")
        testcase.assertEqual(desc["provenance"]["skill"], expected_skill)


class TestClinicalPluginContract(unittest.TestCase):
    def test_plugin_layout_frontmatter_locks_schemas_and_ast(self):
        manifest = (PLUGIN_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("id: clinical", manifest)
        self.assertIn("version: 0.1.0", manifest)
        self.assertIn("toolPrefix: clin", manifest)
        self.assertIn("- agents/Clinical.md", manifest)
        self.assertIn("- skills/prediction-causal-workflow", manifest)
        self.assertIn("- skills/survival-workflow", manifest)

        for req in ("icon.svg", "README.md", "CHANGELOG.md", "agents/Clinical.md"):
            self.assertTrue((PLUGIN_DIR / req).is_file(), f"missing {req}")
        for reserved in ("mcp", "wrappers", "orchestrator"):
            self.assertFalse((PLUGIN_DIR / reserved).exists(), f"reserved dir {reserved} exists")

        agent_text = (PLUGIN_DIR / "agents/Clinical.md").read_text(encoding="utf-8")
        self.assertTrue(agent_text.startswith("---\n"))
        agent_fm = agent_text[4 : agent_text.index("\n---\n", 4)]
        self.assertIn("name: Clinical", agent_fm)
        self.assertIn("- survival-workflow", agent_fm)
        self.assertIn("- prediction-causal-workflow", agent_fm)
        for line in agent_fm.splitlines():
            m = re.match(r"^( +)-\s", line)
            if m:
                self.assertEqual(len(m.group(1)) % 2, 0, f"odd YAML indent in Clinical.md: {line}")

        for plat in ("darwin-arm64", "darwin-x64", "linux-x64"):
            lock = PLUGIN_DIR / f"environments/clinical/locks/{plat}.txt"
            self.assertTrue(lock.is_file(), f"missing lock {plat}")
            ltxt = lock.read_text(encoding="utf-8")
            self.assertIn("@EXPLICIT", ltxt)
            self.assertIn("# spec: resources/plugins/clinical/environments/clinical/environment.yml", ltxt)

        # Verify both skills' SKILL.md contracts, relative links, schemas, and scripts
        skill_specs = [
            (
                SURV_DIR,
                "survival-workflow",
                {"inspect-result.json", "table1-result.json", "km-result.json", "cox-result.json", "time-roc-result.json", "competing-risk-result.json"},
                {"_common.py", "clin_inspect.py", "clin_table1.py", "clin_km.py", "clin_cox.py", "clin_time_roc.py", "clin_competing_risk.py"},
            ),
            (
                PRED_DIR,
                "prediction-causal-workflow",
                {"diagnostic-roc-result.json", "calibration-dca-result.json", "psm-iptw-result.json", "rcs-result.json", "meta-result.json"},
                {"_pred_common.py", "clin_diagnostic_roc.py", "clin_calibration_dca.py", "clin_psm_iptw.py", "clin_rcs.py", "clin_meta.py"},
            ),
        ]

        for sdir, sname, exp_schemas, exp_scripts in skill_specs:
            stext = (sdir / "SKILL.md").read_text(encoding="utf-8")
            self.assertTrue(stext.startswith("---\n"))
            sfm = stext[4 : stext.index("\n---\n", 4)]
            self.assertIn(f"name: {sname}", sfm)
            self.assertIn("  environment: plugin:clinical", sfm)
            self.assertIn("    - Clinical", sfm)
            self.assertIsNone(re.search(r"^\s+toolPrefix:", sfm, re.M))

            for line in sfm.splitlines():
                m = re.match(r"^( +)-\s", line)
                if m:
                    self.assertEqual(len(m.group(1)) % 2, 0, f"odd YAML indent in {sname}/SKILL.md: {line}")

            names = re.findall(r"^    - name:\s+\S+", sfm, re.M)
            outs = re.findall(r"^      output:\s+(\S+)", sfm, re.M)
            self.assertEqual(len(names), len(exp_schemas))
            self.assertEqual(len(outs), len(names))
            for rel_out in outs:
                self.assertTrue((sdir / rel_out).is_file(), f"missing schema {rel_out}")

            for rel_link in re.findall(r"\[[^\]]+\]\(([^)\s]+)\)", stext):
                if not rel_link.startswith(("http://", "https://", "#")):
                    self.assertTrue((sdir / rel_link.split("#")[0]).exists(), f"{sname} missing link {rel_link}")

            actual_schemas = {p.name for p in (sdir / "schemas").glob("*.json")}
            self.assertEqual(actual_schemas, exp_schemas)
            for sp in (sdir / "schemas").glob("*.json"):
                self.assertEqual(json.loads(sp.read_text(encoding="utf-8"))["type"], "object")

            actual_scripts = {p.name for p in (sdir / "scripts").glob("*.py")}
            self.assertEqual(actual_scripts, exp_scripts)
            with tempfile.TemporaryDirectory() as pyc_tmp:
                for py in (sdir / "scripts").glob("*.py"):
                    ast.parse(py.read_text(encoding="utf-8"), filename=str(py))
                    py_compile.compile(str(py), cfile=str(Path(pyc_tmp) / f"{py.stem}.pyc"), doraise=True)

        # Ensure zero __pycache__ directories inside resources/plugins/clinical
        pycaches = list(PLUGIN_DIR.rglob("__pycache__"))
        self.assertEqual(pycaches, [], f"unexpected __pycache__ dirs: {pycaches}")


@unittest.skipUnless(HAS_SCIPY_STACK, "numpy/pandas/scipy/sklearn/matplotlib required for live E2E")
class TestClinicalWorkflowsE2E(unittest.TestCase):
    def test_survival_and_competing_risk_tools(self):
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(2026)
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            n = 100
            samples = [f"PT_{i:03d}" for i in range(1, n + 1)]
            arm = np.array(["Experimental"] * 50 + ["Control"] * 50)
            age = np.round(rng.normal(62.0, 9.0, size=n), 1)
            stage = rng.choice(["Stage_I_II", "Stage_III_IV"], size=n, p=[0.55, 0.45])
            gene_a = np.round(rng.normal(0.0, 1.0, size=n) + np.where(arm == "Control", 0.8, -0.4), 3)
            gene_b = np.round(rng.normal(0.0, 1.0, size=n), 3)

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
                }
            )
            clin_path = tdir / "clinical.tsv"
            clin_df.to_csv(clin_path, sep="\t", index=False)

            mat_df = pd.DataFrame([gene_a, gene_b], index=["GENE_A", "GENE_B"], columns=samples)
            mat_df.index.name = "gene_symbol"
            mat_path = tdir / "expr.tsv"
            mat_df.to_csv(mat_path, sep="\t")

            # 1) clin_inspect
            ins = run_tool(
                SURV_DIR,
                "clin_inspect.py",
                ["--clinical", str(clin_path), "--time_col", "os_months", "--event_col", "os_status", "--matrix", str(mat_path)],
                cwd=tdir,
            )
            validate_output_schema(ins, SURV_DIR / "schemas/inspect-result.json")
            self.assertEqual(ins["n_valid_survival"], n)
            self.assertEqual(ins["inferred_time_unit"], "months")

            # 2) clin_table1
            t1 = run_tool(
                SURV_DIR,
                "clin_table1.py",
                ["--clinical", str(clin_path), "--group_col", "arm", "--variables", "age", "--variables", "stage", "--output_dir", str(tdir / "t1")],
                cwd=tdir,
            )
            validate_output_schema(t1, SURV_DIR / "schemas/table1-result.json")
            assert_artifacts_valid(self, t1, "survival-workflow")

            # 3) clin_km
            km = run_tool(
                SURV_DIR,
                "clin_km.py",
                [
                    "--clinical", str(clin_path), "--time_col", "os_months", "--event_col", "os_status",
                    "--biomarker_col", "GENE_A", "--matrix", str(mat_path),
                    "--landmarks", "12", "--landmarks", "36", "--output_dir", str(tdir / "km"),
                ],
                cwd=tdir,
            )
            validate_output_schema(km, SURV_DIR / "schemas/km-result.json")
            assert_artifacts_valid(self, km, "survival-workflow")
            self.assertLess(km["logrank_pvalue"], 0.05)

            # 4) clin_cox
            cox = run_tool(
                SURV_DIR,
                "clin_cox.py",
                [
                    "--clinical", str(clin_path), "--time_col", "os_months", "--event_col", "os_status",
                    "--covariates", "arm", "--covariates", "stage", "--covariates", "GENE_A",
                    "--matrix", str(mat_path), "--standardize_continuous", "--output_dir", str(tdir / "cox"),
                ],
                cwd=tdir,
            )
            validate_output_schema(cox, SURV_DIR / "schemas/cox-result.json")
            assert_artifacts_valid(self, cox, "survival-workflow")
            self.assertGreater(cox["multivariable"]["c_index"], 0.65)

            # 5) clin_time_roc
            troc = run_tool(
                SURV_DIR,
                "clin_time_roc.py",
                [
                    "--clinical", str(clin_path), "--time_col", "os_months", "--event_col", "os_status",
                    "--feature_cols", "age", "--feature_cols", "GENE_A", "--feature_cols", "GENE_B",
                    "--matrix", str(mat_path), "--landmarks", "12", "--landmarks", "24", "--landmarks", "36",
                    "--output_dir", str(tdir / "troc"),
                ],
                cwd=tdir,
            )
            validate_output_schema(troc, SURV_DIR / "schemas/time-roc-result.json")
            assert_artifacts_valid(self, troc, "survival-workflow")
            self.assertGreater(troc["c_index"], 0.60)

            # 6) clin_competing_risk
            n_cr = 160
            grp_cr = np.array(["HighRisk"] * 80 + ["LowRisk"] * 80)
            age_cr = np.round(rng.normal(64.0, 8.0, size=n_cr), 1)
            # Cause 1 (primary) hazard higher in HighRisk; Cause 2 (competing) substantial in both
            t1_lat = rng.exponential(scale=np.where(grp_cr == "HighRisk", 18.0, 55.0))
            t2_lat = rng.exponential(scale=30.0, size=n_cr)
            c_lat = rng.uniform(24.0, 60.0, size=n_cr)
            obs_cr = np.round(np.maximum(0.2, np.minimum(np.minimum(t1_lat, t2_lat), c_lat)), 2)
            st_cr = np.where((t1_lat <= t2_lat) & (t1_lat <= c_lat), 1, np.where((t2_lat <= t1_lat) & (t2_lat <= c_lat), 2, 0))
            cr_df = pd.DataFrame({"id": [f"CR_{i}" for i in range(n_cr)], "time": obs_cr, "status": st_cr, "risk_group": grp_cr, "age": age_cr})
            cr_path = tdir / "competing.tsv"
            cr_df.to_csv(cr_path, sep="\t", index=False)

            cr = run_tool(
                SURV_DIR,
                "clin_competing_risk.py",
                [
                    "--clinical", str(cr_path), "--time_col", "time", "--status_col", "status",
                    "--primary_code", "1", "--group_col", "risk_group",
                    "--covariates", "risk_group", "--covariates", "age",
                    "--landmarks", "12", "--landmarks", "24", "--output_dir", str(tdir / "cr"),
                ],
                cwd=tdir,
            )
            validate_output_schema(cr, SURV_DIR / "schemas/competing-risk-result.json")
            assert_artifacts_valid(self, cr, "survival-workflow")
            self.assertGreater(cr["n_competing_events"], 15)
            self.assertLessEqual(cr["overall_cif"]["final_cif_primary"] + cr["overall_cif"]["final_cif_competing"], 1.0001)
            self.assertGreater(cr["naive_1_minus_km_overestimation"]["final_overestimation"], 0.02)
            self.assertLess(cr["grays_test"]["pvalue"], 0.05)
            self.assertIsNotNone(cr["fine_gray_model"])
            self.assertGreater(len(cr["landmark_cif"]), 0)

    def test_prediction_causal_and_meta_tools(self):
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(777)
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)

            # 7) clin_diagnostic_roc
            n_d = 180
            y_d = np.array([1] * 70 + [0] * 110)
            strong = rng.normal(loc=np.where(y_d == 1, 1.8, 0.0), scale=1.0)
            weak = rng.normal(loc=np.where(y_d == 1, 0.3, 0.0), scale=1.0)
            inverse = rng.normal(loc=np.where(y_d == 1, -1.5, 0.0), scale=1.0)
            roc_df = pd.DataFrame({"y": y_d, "strong_score": strong, "weak_score": weak, "protective_marker": inverse})
            roc_path = tdir / "diag.tsv"
            roc_df.to_csv(roc_path, sep="\t", index=False)

            droc = run_tool(
                SURV_DIR.parent / "prediction-causal-workflow",
                "clin_diagnostic_roc.py",
                [
                    "--clinical", str(roc_path), "--outcome_col", "y",
                    "--predictors", "strong_score", "--predictors", "weak_score", "--predictors", "protective_marker",
                    "--n_bootstrap", "100", "--output_dir", str(tdir / "droc"),
                ],
                cwd=tdir,
            )
            validate_output_schema(droc, PRED_DIR / "schemas/diagnostic-roc-result.json")
            assert_artifacts_valid(self, droc, "prediction-causal-workflow")
            pmap = {p["predictor"]: p for p in droc["predictors"]}
            self.assertGreater(pmap["strong_score"]["auc"], 0.85)
            self.assertLess(pmap["weak_score"]["auc"], 0.70)
            self.assertTrue(pmap["protective_marker"]["direction_inverted"])
            self.assertGreater(pmap["protective_marker"]["auc"], 0.80)
            dl_map = {d["comparator"]: d for d in droc["delong_comparisons"]}
            self.assertLess(dl_map["weak_score"]["pvalue"], 0.01)
            self.assertLess(dl_map["weak_score"]["delta_auc"], -0.15)

            # 8) clin_calibration_dca (precomputed calibrated probs + OOF CV features)
            n_c = 250
            x1 = rng.normal(0.0, 1.0, size=n_c)
            x2 = rng.normal(0.0, 1.0, size=n_c)
            lin_c = -0.3 + 1.4 * x1 - 0.9 * x2
            p_true = 1.0 / (1.0 + np.exp(-lin_c))
            y_c = rng.binomial(1, p_true)
            cal_df = pd.DataFrame({"y": y_c, "prob_cal": np.round(p_true, 6), "x1": x1, "x2": x2})
            cal_path = tdir / "cal.tsv"
            cal_df.to_csv(cal_path, sep="\t", index=False)

            cal_res = run_tool(
                PRED_DIR,
                "clin_calibration_dca.py",
                ["--clinical", str(cal_path), "--outcome_col", "y", "--prob_col", "prob_cal", "--output_dir", str(tdir / "cal_prob")],
                cwd=tdir,
            )
            validate_output_schema(cal_res, PRED_DIR / "schemas/calibration-dca-result.json")
            assert_artifacts_valid(self, cal_res, "prediction-causal-workflow")
            self.assertGreater(cal_res["calibration"]["calibration_slope"], 0.70)
            self.assertLess(cal_res["calibration"]["calibration_slope"], 1.40)
            self.assertLess(abs(cal_res["calibration"]["citl"]), 0.40)
            self.assertGreater(cal_res["calibration"]["scaled_brier_ipa"], 0.15)
            self.assertIsNotNone(cal_res["dca"]["useful_threshold_range"])
            self.assertGreater(cal_res["dca"]["max_net_benefit"], 0.05)

            cal_oof = run_tool(
                PRED_DIR,
                "clin_calibration_dca.py",
                ["--clinical", str(cal_path), "--outcome_col", "y", "--features", "x1", "--features", "x2", "--output_dir", str(tdir / "cal_oof")],
                cwd=tdir,
            )
            validate_output_schema(cal_oof, PRED_DIR / "schemas/calibration-dca-result.json")
            self.assertEqual(cal_oof["evaluation_mode"], "oof_cv_logistic")

            # 9) clin_psm_iptw (confounded observational cohort, true ATE = +2.0)
            n_p = 240
            age_p = rng.normal(60.0, 10.0, size=n_p)
            sev_p = rng.normal(0.0, 1.0, size=n_p)
            sex_p = rng.choice(["F", "M"], size=n_p)
            logit_t = -0.3 + 0.04 * (age_p - 60.0) + 0.85 * sev_p + 0.3 * (sex_p == "M").astype(float)
            treat_p = rng.binomial(1, 1.0 / (1.0 + np.exp(-logit_t)))
            # True treatment effect = +2.0, while severity negatively impacts outcome (-3.0*sev) so unadjusted is biased down!
            y_p = 10.0 + 2.0 * treat_p - 3.0 * sev_p - 0.08 * (age_p - 60.0) + rng.normal(0.0, 1.0, size=n_p)
            psm_df = pd.DataFrame({"id": [f"P{i}" for i in range(n_p)], "treat": treat_p, "age": age_p, "severity": sev_p, "sex": sex_p, "score": y_p})
            psm_path = tdir / "psm.tsv"
            psm_df.to_csv(psm_path, sep="\t", index=False)

            psm_res = run_tool(
                PRED_DIR,
                "clin_psm_iptw.py",
                [
                    "--clinical", str(psm_path), "--treatment_col", "treat",
                    "--covariates", "age", "--covariates", "severity", "--covariates", "sex",
                    "--outcome_col", "score", "--method", "both", "--caliper_sd", "0.25",
                    "--output_dir", str(tdir / "psm"),
                ],
                cwd=tdir,
            )
            validate_output_schema(psm_res, PRED_DIR / "schemas/psm-iptw-result.json")
            assert_artifacts_valid(self, psm_res, "prediction-causal-workflow")
            sev_bal = [b for b in psm_res["covariate_balance"] if b["covariate"] == "severity"][0]
            self.assertGreater(sev_bal["abs_smd_unadjusted"], 0.40)
            self.assertLess(sev_bal["abs_smd_psm"], 0.15)
            self.assertLess(sev_bal["abs_smd_iptw"], 0.15)
            iptw_est = psm_res["treatment_effects"]["iptw"]["estimate"]
            self.assertAlmostEqual(iptw_est, 2.0, delta=0.6)

            # 10) clin_rcs (U-shaped nonlinear vs linear exposure)
            n_r = 220
            x_r = rng.uniform(0.0, 10.0, size=n_r)
            adj_r = rng.normal(0.0, 1.0, size=n_r)
            y_quad = 0.7 * ((x_r - 5.0) ** 2) + 0.5 * adj_r + rng.normal(0.0, 1.5, size=n_r)
            y_lin = 0.8 * x_r + 0.5 * adj_r + rng.normal(0.0, 1.5, size=n_r)
            rcs_df = pd.DataFrame({"exposure": x_r, "adj": adj_r, "y_quad": y_quad, "y_lin": y_lin})
            rcs_path = tdir / "rcs.tsv"
            rcs_df.to_csv(rcs_path, sep="\t", index=False)

            rcs_nl = run_tool(
                PRED_DIR,
                "clin_rcs.py",
                [
                    "--clinical", str(rcs_path), "--exposure_col", "exposure", "--outcome_col", "y_quad",
                    "--covariates", "adj", "--n_knots", "4", "--ref_value", "5.0", "--output_dir", str(tdir / "rcs_nl"),
                ],
                cwd=tdir,
            )
            validate_output_schema(rcs_nl, PRED_DIR / "schemas/rcs-result.json")
            assert_artifacts_valid(self, rcs_nl, "prediction-causal-workflow")
            self.assertLess(rcs_nl["p_overall"], 0.001)
            self.assertLess(rcs_nl["p_nonlinear"], 0.001)

            rcs_l = run_tool(
                PRED_DIR,
                "clin_rcs.py",
                [
                    "--clinical", str(rcs_path), "--exposure_col", "exposure", "--outcome_col", "y_lin",
                    "--covariates", "adj", "--n_knots", "4", "--output_dir", str(tdir / "rcs_lin"),
                ],
                cwd=tdir,
            )
            validate_output_schema(rcs_l, PRED_DIR / "schemas/rcs-result.json")
            self.assertLess(rcs_l["p_overall"], 0.001)
            self.assertGreater(rcs_l["p_nonlinear"], 0.10)

            # 11) clin_meta (ratio-scale HRs + 2x2 zero-cell OR)
            meta_df = pd.DataFrame(
                {
                    "study": ["Trial_A", "Trial_B", "Trial_C", "Trial_D", "Trial_E"],
                    "hr": [0.65, 0.72, 0.58, 0.88, 0.62],
                    "ci_low": [0.48, 0.55, 0.40, 0.64, 0.45],
                    "ci_high": [0.88, 0.94, 0.84, 1.21, 0.85],
                    "region": ["Asia", "Europe", "Asia", "Europe", "Asia"],
                }
            )
            meta_path = tdir / "meta_hr.tsv"
            meta_df.to_csv(meta_path, sep="\t", index=False)

            m_res = run_tool(
                PRED_DIR,
                "clin_meta.py",
                [
                    "--studies", str(meta_path), "--study_col", "study",
                    "--effect_col", "hr", "--ci_low_col", "ci_low", "--ci_high_col", "ci_high",
                    "--is_ratio_scale", "--subgroup_col", "region", "--output_dir", str(tdir / "meta_hr"),
                ],
                cwd=tdir,
            )
            validate_output_schema(m_res, PRED_DIR / "schemas/meta-result.json")
            assert_artifacts_valid(self, m_res, "prediction-causal-workflow")
            self.assertEqual(m_res["n_studies"], 5)
            # Analytical FE pooled log-HR check
            log_hr = np.log(meta_df["hr"].to_numpy())
            se_hr = (np.log(meta_df["ci_high"].to_numpy()) - np.log(meta_df["ci_low"].to_numpy())) / (2.0 * 1.959964)
            w_fe = 1.0 / (se_hr ** 2)
            exp_fe = math.exp(float(np.sum(w_fe * log_hr) / np.sum(w_fe)))
            self.assertAlmostEqual(m_res["fixed_effect"]["estimate"], exp_fe, places=3)
            self.assertIsNotNone(m_res["prediction_interval_95"])
            self.assertIsNotNone(m_res["eggers_test"])
            self.assertIsNotNone(m_res["subgroups"])
            self.assertEqual(len(m_res["leave_one_out"]), 5)

            # 2x2 contingency mode with zero cell
            tab22 = pd.DataFrame(
                {
                    "study": ["S1", "S2", "S3"],
                    "ai": [12, 0, 18],
                    "bi": [38, 25, 42],
                    "ci": [22, 5, 30],
                    "di": [28, 20, 30],
                }
            )
            tab22_path = tdir / "meta_2x2.tsv"
            tab22.to_csv(tab22_path, sep="\t", index=False)
            m22 = run_tool(
                PRED_DIR,
                "clin_meta.py",
                [
                    "--studies", str(tab22_path), "--study_col", "study",
                    "--ai_col", "ai", "--bi_col", "bi", "--ci_col", "ci", "--di_col", "di",
                    "--measure", "OR", "--output_dir", str(tdir / "meta_22"),
                ],
                cwd=tdir,
            )
            validate_output_schema(m22, PRED_DIR / "schemas/meta-result.json")
            self.assertTrue(any("Haldane-Anscombe" in w for w in m22["warnings"]))


if __name__ == "__main__":
    unittest.main()
