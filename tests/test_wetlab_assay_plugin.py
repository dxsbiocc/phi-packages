#!/usr/bin/env python3
"""Contract, layout, AST, schema, and end-to-end runtime tests for wetlab-assay (0.1.0)."""

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
PLUGIN_DIR = ROOT / "resources" / "plugins" / "wetlab-assay"
SKILL_DIR = PLUGIN_DIR / "skills" / "assay-workflow"

HAS_SCIPY_STACK = all(
    importlib.util.find_spec(pkg) is not None
    for pkg in ("numpy", "pandas", "scipy", "sklearn", "matplotlib")
)
HAS_JSONSCHEMA = importlib.util.find_spec("jsonschema") is not None


def _run_json(script_name: str, args: list[str], cwd: Path) -> dict:
    script = SKILL_DIR / "scripts" / script_name
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
    return json.loads(proc.stdout)


def _validate_schema_and_sidecars(testcase: unittest.TestCase, payload: dict, schema_name: str) -> None:
    schema_path = SKILL_DIR / "schemas" / schema_name
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    if HAS_JSONSCHEMA:
        import jsonschema

        jsonschema.validate(instance=payload, schema=schema)
    for req in schema.get("required", []):
        testcase.assertIn(req, payload)
    testcase.assertGreater(len(payload["artifacts"]), 0)
    for art in payload["artifacts"]:
        p = Path(art)
        sidecar = p.with_name(p.name + ".phi-artifact.json")
        testcase.assertTrue(p.is_file(), f"missing artifact {art}")
        testcase.assertTrue(sidecar.is_file(), f"missing sidecar {sidecar}")
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
        testcase.assertEqual(meta["contractVersion"], "1.0.0")
        testcase.assertEqual(meta["provenance"]["tool"], payload["tool"])
        testcase.assertEqual(meta["provenance"]["skill"], "assay-workflow")


class TestWetlabAssayContract(unittest.TestCase):
    def test_plugin_layout_locks_schemas_and_ast(self) -> None:
        manifest = (PLUGIN_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("id: wetlab-assay", manifest)
        self.assertIn("version: 0.1.0", manifest)
        self.assertIn("toolPrefix: wa", manifest)
        self.assertIn("- agents/WetlabAssay.md", manifest)
        self.assertIn("- skills/assay-workflow", manifest)

        for req in ("icon.svg", "README.md", "CHANGELOG.md", "agents/WetlabAssay.md"):
            self.assertTrue((PLUGIN_DIR / req).is_file(), f"missing {req}")
        for reserved in ("mcp", "wrappers", "orchestrator"):
            self.assertFalse((PLUGIN_DIR / reserved).exists(), f"reserved dir {reserved} exists")

        for plat in ("darwin-arm64", "darwin-x64", "linux-x64"):
            lock = PLUGIN_DIR / f"environments/wa/locks/{plat}.txt"
            self.assertTrue(lock.is_file(), f"missing lock {plat}")
            ltxt = lock.read_text(encoding="utf-8")
            self.assertIn("@EXPLICIT", ltxt)
            self.assertIn("# spec: resources/plugins/wetlab-assay/environments/wa/environment.yml", ltxt)

        stext = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
        sfm = stext[4 : stext.index("\n---\n", 4)]
        self.assertIn("name: assay-workflow", sfm)
        self.assertIn("  environment: plugin:wa", sfm)
        self.assertIn("    - WetlabAssay", sfm)
        self.assertIsNone(re.search(r"^\s+toolPrefix:", sfm, re.M))
        for line in sfm.splitlines():
            m = re.match(r"^( +)-\s", line)
            if m:
                self.assertEqual(len(m.group(1)) % 2, 0, f"odd YAML indent: {line}")

        names = re.findall(r"^    - name:\s+(\S+)", sfm, re.M)
        outs = re.findall(r"^      output:\s+(\S+)", sfm, re.M)
        self.assertEqual(names, ["qpcr_quant", "immunoassay_curve", "enzyme_kinetics", "binding_kinetics"])
        self.assertEqual(len(outs), 4)
        for rel_out in outs:
            self.assertTrue((SKILL_DIR / rel_out).is_file())

        exp_scripts = {
            "_common.py",
            "wa_qpcr_quant.py",
            "wa_immunoassay_curve.py",
            "wa_enzyme_kinetics.py",
            "wa_binding_kinetics.py",
        }
        actual_scripts = {p.name for p in (SKILL_DIR / "scripts").glob("*.py")}
        self.assertEqual(actual_scripts, exp_scripts)
        with tempfile.TemporaryDirectory() as pyc_tmp:
            for py in (SKILL_DIR / "scripts").glob("*.py"):
                ast.parse(py.read_text(encoding="utf-8"), filename=str(py))
                py_compile.compile(str(py), cfile=str(Path(pyc_tmp) / f"{py.stem}.pyc"), doraise=True)

        self.assertEqual(list(PLUGIN_DIR.rglob("__pycache__")), [])


@unittest.skipUnless(HAS_SCIPY_STACK, "numpy/pandas/scipy/sklearn/matplotlib required for live E2E")
class TestWetlabAssayRuntimeE2E(unittest.TestCase):
    def test_1_qpcr_pfaffl_genorm_and_ntc_qc(self) -> None:
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(42)
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # Standard curve: E_IL6 = 1.95 (slope = -1/log10(1.95) = -3.448), ACTB = 2.0, RPLP0 = 1.98
            sc_rows = []
            for gene, eff in (("ACTB", 2.00), ("RPLP0", 1.98), ("IL6", 1.95)):
                slope = -1.0 / math.log10(eff)
                for dil in (1e6, 1e5, 1e4, 1e3, 1e2):
                    for rep in range(2):
                        sc_rows.append(
                            {
                                "gene": gene,
                                "dilution_or_copy": dil,
                                "cq": round(36.0 + slope * math.log10(dil / 100.0) + rng.normal(0.0, 0.03), 3),
                            }
                        )
            sc_path = tdir / "std_curve.tsv"
            pd.DataFrame(sc_rows).to_csv(sc_path, sep="\t", index=False)

            # Biological replicates: 4 Control vs 4 LPS; IL6 shifts by -3.0 Cq (~8-fold up in LPS)
            cq_rows = []
            for grp, shift in (("Control", 0.0), ("LPS", -3.0)):
                for bio in range(1, 5):
                    sname = f"{grp}_{bio}"
                    base_offset = rng.normal(0.0, 0.12)
                    for tech in range(3):
                        cq_rows.append({"sample": sname, "group": grp, "gene": "ACTB", "cq": round(19.0 + base_offset + rng.normal(0.0, 0.05), 3)})
                        cq_rows.append({"sample": sname, "group": grp, "gene": "RPLP0", "cq": round(20.2 + base_offset + rng.normal(0.0, 0.05), 3)})
                        cq_rows.append({"sample": sname, "group": grp, "gene": "IL6", "cq": round(26.5 + shift + base_offset + rng.normal(0.0, 0.08), 3)})
            # NTC wells (Undetermined for ACTB/RPLP0, late 38.5 for IL6)
            cq_rows.append({"sample": "NTC_1", "group": "NTC", "gene": "ACTB", "cq": "Undetermined"})
            cq_rows.append({"sample": "NTC_1", "group": "NTC", "gene": "RPLP0", "cq": "Undetermined"})
            cq_rows.append({"sample": "NTC_1", "group": "NTC", "gene": "IL6", "cq": 38.5})

            cq_path = tdir / "cq.tsv"
            pd.DataFrame(cq_rows).to_csv(cq_path, sep="\t", index=False)

            res = _run_json(
                "wa_qpcr_quant.py",
                [
                    "--cq_table", str(cq_path),
                    "--reference_genes", "ACTB",
                    "--reference_genes", "RPLP0",
                    "--control_group", "Control",
                    "--standard_curve", str(sc_path),
                    "--output_dir", str(tdir / "qpcr_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "qpcr-quant-result.json")
            self.assertEqual(res["n_samples"], 8)
            self.assertEqual(res["n_reference_genes"], 2)
            comp = res["comparisons"][0]
            self.assertGreater(comp["group_mean_fc"], 5.5)
            self.assertLess(comp["p_value"], 0.005)

    def test_2_immunoassay_curve_5pl_and_censoring(self) -> None:
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(123)
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # True 4PL: A=0.05, D=3.20, C=50.0, B=-1.2 -> increasing curve y = D + (A - D)/(1 + (x/C)^1.2)
            def true_sig(c: float) -> float:
                if c <= 0.0:
                    return 0.05
                return 0.05 + (3.20 - 0.05) / (1.0 + (50.0 / c) ** 1.2)

            std_rows = []
            for conc in (0.0, 1.56, 3.125, 6.25, 12.5, 25.0, 50.0, 100.0, 200.0, 400.0):
                for _ in range(3):
                    y = true_sig(conc) * rng.normal(1.0, 0.02) + rng.normal(0.0, 0.002)
                    std_rows.append({"concentration": conc, "signal": round(max(0.01, y), 5)})
            std_path = tdir / "elisa_std.tsv"
            pd.DataFrame(std_rows).to_csv(std_path, sep="\t", index=False)

            unk_rows = [
                {"sample": "U_BlankLike", "signal": 0.051, "dilution_factor": 1.0, "group": "Healthy"},
                {"sample": "U_Mid_Dil2x", "signal": round(true_sig(40.0), 4), "dilution_factor": 2.0, "group": "Patient"},
                {"sample": "U_Saturated", "signal": 3.19, "dilution_factor": 1.0, "group": "Patient"},
            ]
            unk_path = tdir / "elisa_unk.tsv"
            pd.DataFrame(unk_rows).to_csv(unk_path, sep="\t", index=False)

            res = _run_json(
                "wa_immunoassay_curve.py",
                [
                    "--standards", str(std_path),
                    "--unknowns", str(unk_path),
                    "--model", "auto",
                    "--weighting", "1_y2",
                    "--output_dir", str(tdir / "elisa_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "immunoassay-curve-result.json")
            self.assertGreater(res["r2"], 0.99)
            self.assertIsNotNone(res["lloq"])
            self.assertIsNotNone(res["uloq"])
            self.assertLess(res["lloq"], res["uloq"])

            unk_df = pd.read_csv(res["outputs"]["interpolated_unknowns"], sep="\t").set_index("sample")
            self.assertIn(unk_df.loc["U_BlankLike", "range_status"], ("below_LOD", "below_LLOQ"))
            self.assertEqual(unk_df.loc["U_Mid_Dil2x", "range_status"], "in_quantifiable_range")
            # 40 * dilution_factor(2.0) ~= 80.0
            self.assertAlmostEqual(float(unk_df.loc["U_Mid_Dil2x", "reported_conc"]), 80.0, delta=8.0)
            self.assertEqual(unk_df.loc["U_Saturated", "range_status"], "above_ULOQ")

    def test_3_enzyme_kinetics_competitive_and_progress_curves(self) -> None:
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(99)
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # 3a. Competitive inhibition: Vmax=100, Km=10, Kic=5
            rows = []
            for inh in (0.0, 5.0, 15.0):
                km_app = 10.0 * (1.0 + inh / 5.0)
                for s in (2.0, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0):
                    for _ in range(2):
                        v = (100.0 * s) / (km_app + s) * rng.normal(1.0, 0.015)
                        rows.append({"substrate_conc": s, "inhibitor_conc": inh, "velocity": round(v, 4)})
            kin_path = tdir / "comp_inhib.tsv"
            pd.DataFrame(rows).to_csv(kin_path, sep="\t", index=False)

            res_inh = _run_json(
                "wa_enzyme_kinetics.py",
                [
                    "--kinetics_data", str(kin_path),
                    "--enzyme_conc", "0.5",
                    "--ic50", "15.0",
                    "--cheng_prusoff_substrate", "20.0",
                    "--output_dir", str(tdir / "kin_inh"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_inh, "enzyme-kinetics-result.json")
            self.assertTrue(res_inh["has_inhibitor"])
            self.assertEqual(res_inh["best_model"], "Competitive")
            self.assertAlmostEqual(res_inh["kcat"], 200.0, delta=12.0)
            # Ki = IC50 / (1 + [S]/Km) = 15 / (1 + 20/10) = 5.0
            self.assertAlmostEqual(res_inh["cheng_prusoff_ki"], 5.0, delta=0.5)

            # 3b. Progress curves mode: Vmax=50, Km=8
            pc_rows = []
            for s in (2.0, 4.0, 8.0, 16.0, 32.0, 64.0):
                v0 = (50.0 * s) / (8.0 + s)
                for t in (0.0, 10.0, 20.0, 30.0, 40.0, 50.0):
                    # Slight late depletion curvature after 35 s
                    sig = v0 * t * (1.0 - 0.0008 * t)
                    pc_rows.append({"time": t, "substrate_conc": s, "signal": round(sig, 4)})
            pc_path = tdir / "progress.tsv"
            pd.DataFrame(pc_rows).to_csv(pc_path, sep="\t", index=False)

            res_pc = _run_json(
                "wa_enzyme_kinetics.py",
                [
                    "--kinetics_data", str(pc_path),
                    "--from_progress_curves",
                    "--output_dir", str(tdir / "kin_pc"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_pc, "enzyme-kinetics-result.json")
            self.assertEqual(res_pc["mode"], "progress_curves")
            self.assertAlmostEqual(res_pc["best_params"]["Vmax"], 50.0, delta=4.0)
            self.assertAlmostEqual(res_pc["best_params"]["Km"], 8.0, delta=1.5)

    def test_4_binding_kinetics_spr_langmuir(self) -> None:
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(2026)
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # True 1:1 Langmuir: kon = 2e5 M^-1 s^-1, koff = 2e-3 s^-1 -> KD = 10 nM, Rmax = 120 RU, td = 900 s (near steady-state)
            kon_true = 2.0e5
            koff_true = 2.0e-3
            rmax_true = 120.0
            td = 900.0
            times = np.linspace(0.0, 1200.0, 121)
            rows = []
            for conc_nm in (1.25, 2.5, 5.0, 10.0, 20.0, 40.0):
                c_m = conc_nm * 1e-9
                k_obs = kon_true * c_m + koff_true
                r_eq = (kon_true * c_m * rmax_true) / k_obs
                r_td = r_eq * (1.0 - math.exp(-k_obs * td))
                for t in times:
                    if t <= td:
                        r = r_eq * (1.0 - math.exp(-k_obs * t))
                    else:
                        r = r_td * math.exp(-koff_true * (t - td))
                    rows.append(
                        {
                            "time": round(float(t), 2),
                            "conc_nm": conc_nm,
                            "response": round(r + float(rng.normal(0.0, 0.35)), 4),
                        }
                    )
            spr_path = tdir / "spr.tsv"
            pd.DataFrame(rows).to_csv(spr_path, sep="\t", index=False)

            res = _run_json(
                "wa_binding_kinetics.py",
                [
                    "--sensorgrams", str(spr_path),
                    "--dissociation_start", "900",
                    "--conc_unit", "nM",
                    "--output_dir", str(tdir / "spr_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "binding-kinetics-result.json")
            self.assertEqual(res["rmax_mode"], "global")
            self.assertEqual(res["n_concentrations"], 6)
            self.assertAlmostEqual(res["kd_nm"], 10.0, delta=0.8)
            self.assertAlmostEqual(res["rmax"], 120.0, delta=4.0)
            self.assertIsNotNone(res["kd_eq_nm"])
            self.assertAlmostEqual(res["kd_eq_nm"], 10.0, delta=2.5)


if __name__ == "__main__":
    unittest.main()
