#!/usr/bin/env python3
"""Contract, layout, AST, schema, and live synthetic E2E tests for statistical-genetics (0.1.0)."""

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
PLUGIN_DIR = ROOT / "resources" / "plugins" / "statistical-genetics"
SKILL_DIR = PLUGIN_DIR / "skills" / "gwas-workflow"

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
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise AssertionError(f"{script_name} did not emit valid JSON: {exc}\nSTDOUT:\n{proc.stdout}") from exc


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
        testcase.assertEqual(meta["provenance"]["skill"], "gwas-workflow")


class TestStatisticalGeneticsContract(unittest.TestCase):
    def test_plugin_layout_locks_schemas_and_ast(self) -> None:
        manifest = (PLUGIN_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("schemaVersion: 1", manifest)
        self.assertIn("id: statistical-genetics", manifest)
        self.assertIn("type: plugin", manifest)
        self.assertIn("version: 0.1.0", manifest)
        self.assertIn("title: Statistical Genetics", manifest)
        self.assertIn("toolPrefix: gwas", manifest)
        self.assertIn("- agents/StatisticalGenetics.md", manifest)
        self.assertIn("- skills/gwas-workflow", manifest)
        self.assertIn("spec: environments/gwas/environment.yml", manifest)

        m_sum = re.search(r"^summary:\s*>-\n((?:  .*\n)+)", manifest, re.M)
        self.assertIsNotNone(m_sum)
        summary_text = " ".join(line.strip() for line in m_sum.group(1).splitlines())
        self.assertLessEqual(len(summary_text), 260, f"summary too long: {len(summary_text)}")

        for req in ("icon.svg", "README.md", "CHANGELOG.md", "agents/StatisticalGenetics.md"):
            self.assertTrue((PLUGIN_DIR / req).is_file(), f"missing {req}")
        for reserved in ("mcp", "wrappers", "orchestrator"):
            self.assertFalse((PLUGIN_DIR / reserved).exists(), f"reserved dir {reserved} exists")

        readme_txt = (PLUGIN_DIR / "README.md").read_text(encoding="utf-8")
        for wname in ("plink2", "regenie", "gcta", "beagle5", "shapeit5", "mendelian-randomization"):
            self.assertIn(wname, readme_txt)

        agent_txt = (PLUGIN_DIR / "agents/StatisticalGenetics.md").read_text(encoding="utf-8")
        for wtool in (
            "wrapper_search",
            "wrapper_inspect",
            "wrapper_plan_run",
            "wrapper_run",
            "wrapper_status",
            "wrapper_wait",
            "wrapper_cancel",
        ):
            self.assertIn(f"- {wtool}", agent_txt)

        for plat in ("darwin-arm64", "darwin-x64", "linux-x64"):
            lock = PLUGIN_DIR / f"environments/gwas/locks/{plat}.txt"
            self.assertTrue(lock.is_file(), f"missing lock {plat}")
            ltxt = lock.read_text(encoding="utf-8")
            self.assertIn("@EXPLICIT", ltxt)
            self.assertIn("# spec: resources/plugins/statistical-genetics/environments/gwas/environment.yml", ltxt)

        stext = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
        sfm = stext[4 : stext.index("\n---\n", 4)]
        sbody = stext[stext.index("\n---\n", 4) + 5 :]
        self.assertIn("name: gwas-workflow", sfm)
        self.assertIn("  environment: plugin:gwas", sfm)
        self.assertIn("    - StatisticalGenetics", sfm)
        self.assertIsNone(re.search(r"^\s*toolPrefix:", sfm, re.M))
        self.assertIn("[references/analysis-playbook.md](references/analysis-playbook.md)", sbody)
        self.assertTrue((SKILL_DIR / "references" / "analysis-playbook.md").is_file())

        for line in sfm.splitlines():
            m = re.match(r"^( +)-\s", line)
            if m:
                self.assertEqual(len(m.group(1)) % 2, 0, f"odd YAML indent: {line}")

        names = re.findall(r"^    - name:\s+(\S+)", sfm, re.M)
        outs = re.findall(r"^      output:\s+(\S+)", sfm, re.M)
        self.assertEqual(
            names,
            [
                "sumstats_harmonize",
                "ldsc_interpret",
                "finemap",
                "meta_qc",
                "prs_eval",
            ],
        )
        self.assertEqual(len(outs), 5)
        for rel_out in outs:
            schema_file = SKILL_DIR / rel_out
            self.assertTrue(schema_file.is_file(), f"missing schema {rel_out}")
            json.loads(schema_file.read_text(encoding="utf-8"))

        exp_scripts = {
            "_common.py",
            "gwas_sumstats_harmonize.py",
            "gwas_ldsc_interpret.py",
            "gwas_finemap.py",
            "gwas_meta_qc.py",
            "gwas_prs_eval.py",
        }
        actual_scripts = {p.name for p in (SKILL_DIR / "scripts").glob("*.py")}
        self.assertEqual(actual_scripts, exp_scripts)
        with tempfile.TemporaryDirectory() as pyc_tmp:
            for py in (SKILL_DIR / "scripts").glob("*.py"):
                src = py.read_text(encoding="utf-8")
                ast.parse(src, filename=str(py))
                py_compile.compile(str(py), cfile=str(Path(pyc_tmp) / f"{py.stem}.pyc"), doraise=True)
                self.assertIn('os.environ.setdefault("MPLBACKEND", "Agg")', src)

        self.assertEqual(list(PLUGIN_DIR.rglob("__pycache__")), [])

    def test_error_json_on_missing_input(self) -> None:
        if not HAS_SCIPY_STACK:
            self.skipTest("scipy stack required to import script")
        with tempfile.TemporaryDirectory() as tmp:
            proc = subprocess.run(
                [sys.executable, "-B", str(SKILL_DIR / "scripts/gwas_sumstats_harmonize.py"), "--sumstats", "missing.tsv", "--output_dir", tmp],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 1)
            err = json.loads(proc.stdout)
            self.assertIn("error", err)


@unittest.skipUnless(HAS_SCIPY_STACK, "numpy/pandas/scipy/sklearn/matplotlib required for live E2E")
class TestStatisticalGeneticsRuntimeE2E(unittest.TestCase):
    def test_1_sumstats_harmonize_regenie_and_or(self) -> None:
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(101)
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # REGENIE-like input with LOG10P, ALLELE0/ALLELE1, A1FREQ, INFO, plus bad variants
            rows = []
            for i in range(120):
                ch = 1 + (i % 4)
                pos = 1_000_000 + (i // 4) * 50_000
                z = float(rng.normal(0.0, 1.0))
                if i == 10:
                    ch, pos, z = 1, 2_000_000, 7.2
                elif i == 11:
                    ch, pos, z = 1, 2_100_000, 6.1  # within 500 kb of SNP 10 -> clumped
                elif i == 50:
                    ch, pos, z = 2, 5_000_000, -6.5  # independent 2nd lead locus
                se = 0.02
                beta = z * se
                log10p = -math.log10(max(2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2.0)))), 1e-300))
                rows.append(
                    {
                        "CHROM": ch,
                        "GENPOS": pos,
                        "ID": f"rs{1000 + i}",
                        "ALLELE0": "A",
                        "ALLELE1": "G",
                        "A1FREQ": 0.25,
                        "INFO": 0.95,
                        "BETA": beta,
                        "SE": se,
                        "LOG10P": log10p,
                    }
                )
            # Add low MAF, low INFO, invalid allele, and SE <= 0 rows
            rows.append({"CHROM": 3, "GENPOS": 900_000, "ID": "rs_lowmaf", "ALLELE0": "C", "ALLELE1": "T", "A1FREQ": 0.002, "INFO": 0.95, "BETA": 0.05, "SE": 0.02, "LOG10P": 1.0})
            rows.append({"CHROM": 3, "GENPOS": 950_000, "ID": "rs_lowinfo", "ALLELE0": "C", "ALLELE1": "T", "A1FREQ": 0.20, "INFO": 0.45, "BETA": 0.05, "SE": 0.02, "LOG10P": 1.0})
            rows.append({"CHROM": 3, "GENPOS": 960_000, "ID": "rs_badallele", "ALLELE0": "N", "ALLELE1": "T", "A1FREQ": 0.20, "INFO": 0.95, "BETA": 0.05, "SE": 0.02, "LOG10P": 1.0})
            rows.append({"CHROM": 3, "GENPOS": 970_000, "ID": "rs_badse", "ALLELE0": "C", "ALLELE1": "T", "A1FREQ": 0.20, "INFO": 0.95, "BETA": 0.05, "SE": 0.0, "LOG10P": 1.0})

            reg_path = tdir / "regenie.tsv"
            pd.DataFrame(rows).to_csv(reg_path, sep="\t", index=False)

            out_dir = tdir / "out_regenie"
            res = _run_json(
                "gwas_sumstats_harmonize.py",
                [
                    "--sumstats", str(reg_path),
                    "--output_dir", str(out_dir),
                    "--n_cases", "5000",
                    "--n_controls", "15000",
                    "--min_maf", "0.01",
                    "--min_info", "0.80",
                    "--lead_window_kb", "500",
                ],
                tdir,
            )
            _validate_schema_and_sidecars(self, res, "sumstats-harmonize-result.json")
            self.assertEqual(res["n_variants_input"], 124)
            self.assertEqual(res["n_variants_harmonized"], 120)
            self.assertEqual(res["n_dropped_qc"], 4)
            self.assertEqual(res["n_lead_loci"], 2)
            self.assertAlmostEqual(res["n_eff"], 15000.0, places=1)
            self.assertIsNotNone(res["lambda_1000"])
            self.assertTrue((out_dir / "harmonized_sumstats.tsv").is_file())
            self.assertTrue((out_dir / "lead_loci.tsv").is_file())
            self.assertTrue((out_dir / "gwas_manhattan_qq.png").is_file())

            # OR-scale input test
            or_path = tdir / "or_gwas.csv"
            pd.DataFrame(
                [
                    {"rsid": "rs1", "chr": "1", "pos": 1000, "effect_allele": "A", "other_allele": "G", "eaf": 0.30, "OR": math.exp(0.25), "SE": 0.03, "P": 1e-16, "N": 10000},
                    {"rsid": "rs2", "chr": "2", "pos": 2000, "effect_allele": "C", "other_allele": "T", "eaf": 0.40, "OR": math.exp(-0.04), "SE": 0.03, "P": 0.18, "N": 10000},
                ]
            ).to_csv(or_path, index=False)
            res_or = _run_json("gwas_sumstats_harmonize.py", ["--sumstats", str(or_path), "--output_dir", str(tdir / "out_or")], tdir)
            _validate_schema_and_sidecars(self, res_or, "sumstats-harmonize-result.json")
            self.assertEqual(res_or["effect_scale_input"], "or_converted_to_log_odds")
            self.assertEqual(res_or["n_lead_loci"], 1)

    def test_2_ldsc_interpret_univariate_bivariate_and_proxy(self) -> None:
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(202)
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            m = 4000
            n1, n2 = 3_000.0, 2_800.0
            h2_1_true, h2_2_true, rg_true = 0.32, 0.26, 0.65
            gcov_true = rg_true * math.sqrt(h2_1_true * h2_2_true)

            l2 = rng.uniform(1.0, 20.0, size=m)
            # Generate correlated Z1, Z2 matching E[Z1^2] = 1 + N1*h2_1*l2/M and E[Z1*Z2] = sqrt(N1*N2)*gcov*l2/M
            z1 = np.empty(m, dtype=float)
            z2 = np.empty(m, dtype=float)
            for j in range(m):
                v1 = 1.0 + (n1 * h2_1_true * l2[j]) / m
                v2 = 1.0 + (n2 * h2_2_true * l2[j]) / m
                c12 = (math.sqrt(n1 * n2) * gcov_true * l2[j]) / m
                draw = rng.multivariate_normal([0.0, 0.0], [[v1, c12], [c12, v2]])
                z1[j], z2[j] = draw[0], draw[1]

            ea1 = np.array(["A"] * m)
            nea1 = np.array(["G"] * m)
            ea2 = ea1.copy()
            nea2 = nea1.copy()
            # Swap alleles and negate Z2 on every 5th SNP to test allele alignment
            swap_mask = (np.arange(m) % 5) == 0
            ea2[swap_mask] = "G"
            nea2[swap_mask] = "A"
            z2[swap_mask] = -z2[swap_mask]

            t1_df = pd.DataFrame(
                {
                    "SNP": [f"rs{j}" for j in range(m)],
                    "CHR": [(j // 100) + 1 for j in range(m)],
                    "POS": [(j % 100) * 50_000 + 100_000 for j in range(m)],
                    "EA": ea1,
                    "NEA": nea1,
                    "Z": z1,
                    "N": n1,
                    "l2": l2,
                }
            )
            t2_df = pd.DataFrame(
                {
                    "SNP": [f"rs{j}" for j in range(m)],
                    "CHR": [(j // 100) + 1 for j in range(m)],
                    "POS": [(j % 100) * 50_000 + 100_000 for j in range(m)],
                    "EA": ea2,
                    "NEA": nea2,
                    "Z": z2,
                    "N": n2,
                }
            )
            p1 = tdir / "trait1.tsv"
            p2 = tdir / "trait2.tsv"
            t1_df.to_csv(p1, sep="\t", index=False)
            t2_df.to_csv(p2, sep="\t", index=False)

            out_biv = tdir / "out_biv"
            res = _run_json(
                "gwas_ldsc_interpret.py",
                [
                    "--trait1_sumstats", str(p1),
                    "--trait2_sumstats", str(p2),
                    "--sample_prev", "0.30",
                    "--pop_prev", "0.05",
                    "--n_blocks", "30",
                    "--output_dir", str(out_biv),
                ],
                tdir,
            )
            _validate_schema_and_sidecars(self, res, "ldsc-interpret-result.json")
            self.assertEqual(res["ld_score_source"], "embedded_l2_column")
            self.assertTrue(res["bivariate_executed"])
            self.assertGreater(res["h2_observed"], 0.15)
            self.assertLess(res["h2_observed"], 0.50)
            self.assertIsNotNone(res["h2_liability"])
            self.assertGreater(res["h2_liability"], res["h2_observed"])
            self.assertGreater(res["ldsc_intercept"], 0.70)
            self.assertLess(res["ldsc_intercept"], 1.35)
            self.assertIsNotNone(res["rg"])
            self.assertGreater(res["rg"], 0.35)
            self.assertLess(res["rg"], 0.95)
            self.assertLess(res["rg_p"], 1e-3)

            # Proxy l2 fallback when l2 is omitted
            p1_nol2 = tdir / "trait1_nol2.tsv"
            t1_df.drop(columns=["l2"]).to_csv(p1_nol2, sep="\t", index=False)
            res_proxy = _run_json(
                "gwas_ldsc_interpret.py",
                ["--trait1_sumstats", str(p1_nol2), "--n_blocks", "20", "--output_dir", str(tdir / "out_proxy")],
                tdir,
            )
            _validate_schema_and_sidecars(self, res_proxy, "ldsc-interpret-result.json")
            self.assertEqual(res_proxy["ld_score_source"], "local_window_proxy")
            self.assertTrue(any("local physical-window" in w for w in res_proxy["warnings"]))

    def test_3_finemap_susie_two_causal_and_abf_fallback(self) -> None:
        import numpy as np
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            m = 60
            n = 25_000.0
            # Build 2 independent LD blocks of 30 SNPs each (AR(1) rho=0.85)
            r_mat = np.eye(m, dtype=float)
            for block_start in (0, 30):
                for i in range(block_start, block_start + 30):
                    for j in range(block_start, block_start + 30):
                        r_mat[i, j] = 0.85 ** abs(i - j)

            causal_1, causal_2 = 10, 45
            beta_true = np.zeros(m, dtype=float)
            beta_true[causal_1] = 0.09
            beta_true[causal_2] = -0.08
            z_scores = math.sqrt(n) * (r_mat @ beta_true)
            p_vals = np.maximum(2.0 * np.exp(-0.5 * z_scores * z_scores) / np.sqrt(2.0 * np.pi * np.maximum(z_scores * z_scores, 1.0)), 1e-300)

            locus_df = pd.DataFrame(
                {
                    "SNP": [f"rs_locus_{i}" for i in range(m)],
                    "CHR": ["6"] * m,
                    "POS": [30_000_000 + i * 5_000 for i in range(m)],
                    "Z": z_scores,
                    "P": p_vals,
                    "N": [n] * m,
                }
            )
            locus_path = tdir / "locus.tsv"
            ld_path = tdir / "ld.txt"
            locus_df.to_csv(locus_path, sep="\t", index=False)
            np.savetxt(ld_path, r_mat, fmt="%.6f")

            out_susie = tdir / "out_susie"
            res = _run_json(
                "gwas_finemap.py",
                [
                    "--locus_sumstats", str(locus_path),
                    "--ld_matrix", str(ld_path),
                    "--max_causal", "5",
                    "--coverage", "0.95",
                    "--min_abs_corr", "0.50",
                    "--output_dir", str(out_susie),
                ],
                tdir,
            )
            _validate_schema_and_sidecars(self, res, "finemap-result.json")
            self.assertEqual(res["method"], "SuSiE-RSS")
            self.assertEqual(res["n_credible_sets"], 2, f"Expected 2 credible sets, got {res['credible_sets']}")

            pip_df = pd.read_csv(out_susie / "finemap_pip.tsv", sep="\t")
            row_c1 = pip_df[pip_df["SNP"] == f"rs_locus_{causal_1}"].iloc[0]
            row_c2 = pip_df[pip_df["SNP"] == f"rs_locus_{causal_2}"].iloc[0]
            self.assertGreater(float(row_c1["pip"]), 0.80)
            self.assertGreater(float(row_c2["pip"]), 0.80)
            self.assertIn(int(row_c1["cs_id"]), (1, 2))
            self.assertIn(int(row_c2["cs_id"]), (1, 2))
            self.assertNotEqual(int(row_c1["cs_id"]), int(row_c2["cs_id"]))

            # Fallback without ld_matrix -> Wakefield ABF L=1
            res_abf = _run_json(
                "gwas_finemap.py",
                ["--locus_sumstats", str(locus_path), "--output_dir", str(tdir / "out_abf")],
                tdir,
            )
            _validate_schema_and_sidecars(self, res_abf, "finemap-result.json")
            self.assertEqual(res_abf["method"], "Wakefield_ABF_L1")
            self.assertEqual(res_abf["n_credible_sets"], 1)
            self.assertTrue(any("Wakefield ABF" in w for w in res_abf["warnings"]))

    def test_4_meta_qc_allele_flip_eaf_outlier_and_heterogeneity(self) -> None:
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            c1 = pd.DataFrame(
                [
                    {"SNP": "rs_concord", "CHR": 1, "POS": 1000, "EA": "A", "NEA": "G", "EAF": 0.30, "BETA": 0.14, "SE": 0.02, "N": 10000},
                    {"SNP": "rs_flip", "CHR": 1, "POS": 2000, "EA": "C", "NEA": "T", "EAF": 0.25, "BETA": 0.12, "SE": 0.02, "N": 10000},
                    {"SNP": "rs_strand", "CHR": 1, "POS": 3000, "EA": "A", "NEA": "C", "EAF": 0.20, "BETA": 0.10, "SE": 0.02, "N": 10000},
                    {"SNP": "rs_incompat", "CHR": 1, "POS": 4000, "EA": "A", "NEA": "G", "EAF": 0.30, "BETA": 0.05, "SE": 0.02, "N": 10000},
                    {"SNP": "rs_pal_amb", "CHR": 1, "POS": 5000, "EA": "A", "NEA": "T", "EAF": 0.49, "BETA": 0.06, "SE": 0.02, "N": 10000},
                    {"SNP": "rs_eaf_out", "CHR": 1, "POS": 6000, "EA": "G", "NEA": "T", "EAF": 0.12, "BETA": 0.08, "SE": 0.02, "N": 10000},
                    {"SNP": "rs_het", "CHR": 1, "POS": 7000, "EA": "A", "NEA": "C", "EAF": 0.35, "BETA": 0.20, "SE": 0.02, "N": 10000},
                ]
            )
            c2 = pd.DataFrame(
                [
                    {"SNP": "rs_concord", "CHR": 1, "POS": 1000, "EA": "A", "NEA": "G", "EAF": 0.32, "BETA": 0.15, "SE": 0.02, "N": 12000},
                    # Swapped alleles + negated BETA + flipped EAF
                    {"SNP": "rs_flip", "CHR": 1, "POS": 2000, "EA": "T", "NEA": "C", "EAF": 0.74, "BETA": -0.11, "SE": 0.02, "N": 12000},
                    # Strand flipped A/C -> T/G
                    {"SNP": "rs_strand", "CHR": 1, "POS": 3000, "EA": "T", "NEA": "G", "EAF": 0.21, "BETA": 0.11, "SE": 0.02, "N": 12000},
                    # Incompatible A/G vs A/C -> dropped
                    {"SNP": "rs_incompat", "CHR": 1, "POS": 4000, "EA": "A", "NEA": "C", "EAF": 0.30, "BETA": 0.05, "SE": 0.02, "N": 12000},
                    # Ambiguous palindromic near 0.5 -> dropped
                    {"SNP": "rs_pal_amb", "CHR": 1, "POS": 5000, "EA": "A", "NEA": "T", "EAF": 0.50, "BETA": 0.06, "SE": 0.02, "N": 12000},
                    # EAF outlier (0.12 vs 0.55)
                    {"SNP": "rs_eaf_out", "CHR": 1, "POS": 6000, "EA": "G", "NEA": "T", "EAF": 0.55, "BETA": 0.07, "SE": 0.02, "N": 12000},
                    # High heterogeneity (+0.20 vs -0.15)
                    {"SNP": "rs_het", "CHR": 1, "POS": 7000, "EA": "A", "NEA": "C", "EAF": 0.36, "BETA": -0.15, "SE": 0.02, "N": 12000},
                ]
            )
            p1, p2 = tdir / "cohort1.tsv", tdir / "cohort2.tsv"
            c1.to_csv(p1, sep="\t", index=False)
            c2.to_csv(p2, sep="\t", index=False)

            out_meta = tdir / "out_meta"
            res = _run_json(
                "gwas_meta_qc.py",
                [
                    "--cohorts", str(p1), str(p2),
                    "--cohort_names", "UKB", "BioBankJP",
                    "--eaf_diff_max", "0.20",
                    "--output_dir", str(out_meta),
                ],
                tdir,
            )
            _validate_schema_and_sidecars(self, res, "meta-qc-result.json")
            self.assertEqual(res["n_cohorts"], 2)
            self.assertEqual(res["n_alleles_flipped"], 1)
            self.assertEqual(res["n_incompatible_dropped"], 1)
            self.assertEqual(res["n_palindromic_ambiguous_dropped"], 1)
            self.assertEqual(res["n_eaf_outliers"], 1)
            self.assertGreaterEqual(res["n_high_heterogeneity"], 1)

            mdf = pd.read_csv(out_meta / "gwas_meta_sumstats.tsv", sep="\t").set_index("SNP")
            self.assertEqual(mdf.loc["rs_flip", "Direction"], "++")
            self.assertGreater(float(mdf.loc["rs_het", "I2"]), 90.0)

    def test_5_prs_eval_binary_and_survival(self) -> None:
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(505)
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # 1) Binary case-control cohort
            n_bin = 600
            age = rng.normal(55.0, 9.0, size=n_bin)
            sex = rng.integers(0, 2, size=n_bin)
            pc1 = rng.normal(0.0, 1.0, size=n_bin)
            prs = rng.normal(0.0, 1.0, size=n_bin)
            logit = -1.2 + 0.04 * (age - 55.0) + 0.5 * sex + 0.85 * prs
            prob = 1.0 / (1.0 + np.exp(-logit))
            case = (rng.uniform(0.0, 1.0, size=n_bin) < prob).astype(int)

            bin_path = tdir / "prs_binary.tsv"
            pd.DataFrame({"prs_raw": prs, "cad": case, "age": age, "sex": sex, "PC1": pc1}).to_csv(bin_path, sep="\t", index=False)
            out_bin = tdir / "out_prs_bin"
            res_bin = _run_json(
                "gwas_prs_eval.py",
                [
                    "--cohort", str(bin_path),
                    "--prs_col", "prs_raw",
                    "--outcome_col", "cad",
                    "--covariates", "age", "sex", "PC1",
                    "--pop_prev", "0.08",
                    "--output_dir", str(out_bin),
                ],
                tdir,
            )
            _validate_schema_and_sidecars(self, res_bin, "prs-eval-result.json")
            self.assertEqual(res_bin["outcome_mode"], "binary")
            self.assertEqual(res_bin["prs_effect_measure"], "OR")
            self.assertGreater(res_bin["prs_effect_per_sd"], 1.5)
            self.assertGreater(res_bin["model_comparison"]["delta_r2"], 0.02)
            self.assertGreater(res_bin["model_comparison"]["delta_discrimination"], 0.02)
            self.assertIsNotNone(res_bin["model_comparison"]["r2_liability_delta"])
            self.assertLess(res_bin["model_comparison"]["lrt_p_value"], 1e-4)

            qmap = {r["stratum"]: r for r in res_bin["quantile_risk"]}
            self.assertEqual(qmap["Middle 40-60%"]["estimate"], 1.0)
            self.assertGreater(qmap["Top 10%"]["estimate"], 1.5)
            self.assertLess(qmap["Bottom 10%"]["estimate"], 1.0)

            # 2) Survival Cox PH cohort
            n_surv = 500
            age_s = rng.normal(60.0, 8.0, size=n_surv)
            prs_s = rng.normal(0.0, 1.0, size=n_surv)
            lin_haz = 0.03 * (age_s - 60.0) + 0.65 * prs_s
            event_t = rng.exponential(scale=10.0 * np.exp(-lin_haz))
            censor_t = rng.uniform(3.0, 12.0, size=n_surv)
            obs_t = np.maximum(np.minimum(event_t, censor_t), 0.05)
            event_ind = (event_t <= censor_t).astype(int)

            surv_path = tdir / "prs_survival.tsv"
            pd.DataFrame({"prs_score": prs_s, "event": event_ind, "time_yrs": obs_t, "age": age_s}).to_csv(surv_path, sep="\t", index=False)
            out_surv = tdir / "out_prs_surv"
            res_surv = _run_json(
                "gwas_prs_eval.py",
                [
                    "--cohort", str(surv_path),
                    "--prs_col", "prs_score",
                    "--outcome_col", "event",
                    "--time_col", "time_yrs",
                    "--covariates", "age",
                    "--output_dir", str(out_surv),
                ],
                tdir,
            )
            _validate_schema_and_sidecars(self, res_surv, "prs-eval-result.json")
            self.assertEqual(res_surv["outcome_mode"], "survival")
            self.assertEqual(res_surv["prs_effect_measure"], "HR")
            self.assertEqual(res_surv["model_comparison"]["discrimination_metric"], "Harrell_C_index")
            self.assertGreater(res_surv["prs_effect_per_sd"], 1.3)
            self.assertGreater(res_surv["model_comparison"]["delta_discrimination"], 0.02)
            self.assertLess(res_surv["model_comparison"]["lrt_p_value"], 1e-4)


if __name__ == "__main__":
    unittest.main()
