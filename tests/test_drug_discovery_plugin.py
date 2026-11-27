#!/usr/bin/env python3
"""Contract, layout, AST, schema, pure-Python SMILES graph, and live E2E tests for drug-discovery (0.1.0)."""

from __future__ import annotations

import ast
import csv
import importlib.util
import json
import py_compile
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLUGIN_DIR = ROOT / "resources" / "plugins" / "drug-discovery"
SKILL_DIR = PLUGIN_DIR / "skills" / "chem-workflow"
SCRIPTS_DIR = SKILL_DIR / "scripts"

HAS_JSONSCHEMA = importlib.util.find_spec("jsonschema") is not None


def _run_json(script_name: str, args: list[str], cwd: Path, expect_fail: bool = False) -> dict:
    script = SCRIPTS_DIR / script_name
    proc = subprocess.run(
        [sys.executable, "-B", str(script), *args],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )
    if expect_fail:
        if proc.returncode == 0:
            raise AssertionError(f"{script_name} unexpectedly succeeded:\nSTDOUT:\n{proc.stdout}")
        return json.loads(proc.stdout)
    if proc.returncode != 0:
        raise AssertionError(
            f"{script_name} failed (rc={proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


def _read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _validate_schema_and_sidecars(testcase: unittest.TestCase, payload: dict, schema_name: str) -> None:
    schema_path = SKILL_DIR / "schemas" / schema_name
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    if HAS_JSONSCHEMA:
        import jsonschema

        jsonschema.validate(instance=payload, schema=schema)
    for req in schema.get("required", []):
        testcase.assertIn(req, payload, f"missing required key '{req}' in {schema_name}")
    testcase.assertIn(payload["backend"], ("rdkit", "pure_python_smiles"))
    testcase.assertGreater(len(payload["artifacts"]), 0)
    for art in payload["artifacts"]:
        p = Path(art)
        sidecar = p.with_name(p.name + ".phi-artifact.json")
        testcase.assertTrue(p.is_file(), f"missing artifact {art}")
        testcase.assertTrue(sidecar.is_file(), f"missing sidecar {sidecar}")
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
        testcase.assertEqual(meta["contractVersion"], "1.0.0")
        testcase.assertEqual(meta["provenance"]["tool"], payload["tool"])
        testcase.assertEqual(meta["provenance"]["skill"], "chem-workflow")
        if p.suffix == ".png":
            raw = p.read_bytes()
            testcase.assertTrue(raw.startswith(b"\x89PNG\r\n\x1a\n"))
            testcase.assertGreater(meta["figure"]["widthPx"], 0)
            testcase.assertGreater(meta["figure"]["heightPx"], 0)


class TestDrugDiscoveryContract(unittest.TestCase):
    def test_plugin_layout_locks_schemas_and_ast(self) -> None:
        manifest = (PLUGIN_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("schemaVersion: 1", manifest)
        self.assertIn("id: drug-discovery", manifest)
        self.assertIn("type: plugin", manifest)
        self.assertIn("version: 0.1.0", manifest)
        self.assertIn("toolPrefix: chem", manifest)
        self.assertIn("- agents/DrugDiscovery.md", manifest)
        self.assertIn("- skills/chem-workflow", manifest)

        # Verify summary <= 260 chars
        m_sum = re.search(r"^summary:\s*>-\n((?:  .+\n)+)", manifest, re.M)
        self.assertIsNotNone(m_sum)
        summary_text = " ".join(ln.strip() for ln in m_sum.group(1).splitlines())
        self.assertLessEqual(len(summary_text), 260, f"summary too long ({len(summary_text)} chars)")

        for req in ("icon.svg", "README.md", "CHANGELOG.md", "agents/DrugDiscovery.md"):
            self.assertTrue((PLUGIN_DIR / req).is_file(), f"missing {req}")
        for reserved in ("mcp", "wrappers", "orchestrator"):
            self.assertFalse((PLUGIN_DIR / reserved).exists(), f"reserved dir {reserved} exists")

        readme = (PLUGIN_DIR / "README.md").read_text(encoding="utf-8")
        for ref in ("local/boltz2", "local/chai1", "local/protenix", "local/alphafold3", "pharmacogenomics"):
            self.assertIn(ref, readme)

        for plat in ("darwin-arm64", "darwin-x64", "linux-x64"):
            lock = PLUGIN_DIR / f"environments/chem/locks/{plat}.txt"
            self.assertTrue(lock.is_file(), f"missing lock {plat}")
            ltxt = lock.read_text(encoding="utf-8")
            self.assertIn("@EXPLICIT", ltxt)
            self.assertIn("# spec: resources/plugins/drug-discovery/environments/chem/environment.yml", ltxt)

        stext = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
        sfm = stext[4 : stext.index("\n---\n", 4)]
        self.assertIn("name: chem-workflow", sfm)
        self.assertIn("  environment: plugin:chem", sfm)
        self.assertIn("    - DrugDiscovery", sfm)
        self.assertIsNone(re.search(r"^\s+toolPrefix:", sfm, re.M))
        for line in sfm.splitlines():
            m = re.match(r"^( +)-\s", line)
            if m:
                self.assertEqual(len(m.group(1)) % 2, 0, f"odd YAML indent: {line}")

        names = re.findall(r"^    - name:\s+(\S+)", sfm, re.M)
        outs = re.findall(r"^      output:\s+(\S+)", sfm, re.M)
        self.assertEqual(names, ["standardize", "medchem_triage", "sar_cliff", "admet_profile", "docking_eval"])
        self.assertEqual(len(outs), 5)
        for rel_out in outs:
            self.assertTrue((SKILL_DIR / rel_out).is_file(), f"missing schema {rel_out}")
            json.loads((SKILL_DIR / rel_out).read_text(encoding="utf-8"))

        self.assertIn("[references/analysis-playbook.md](references/analysis-playbook.md)", stext)
        self.assertTrue((SKILL_DIR / "references" / "analysis-playbook.md").is_file())

        exp_scripts = {
            "_common.py",
            "chem_standardize.py",
            "chem_medchem_triage.py",
            "chem_sar_cliff.py",
            "chem_admet_profile.py",
            "chem_docking_eval.py",
        }
        actual_scripts = {p.name for p in SCRIPTS_DIR.glob("*.py")}
        self.assertEqual(actual_scripts, exp_scripts)
        with tempfile.TemporaryDirectory() as pyc_tmp:
            for py in SCRIPTS_DIR.glob("*.py"):
                ast.parse(py.read_text(encoding="utf-8"), filename=str(py))
                py_compile.compile(str(py), cfile=str(Path(pyc_tmp) / f"{py.stem}.pyc"), doraise=True)

        self.assertEqual(list(PLUGIN_DIR.rglob("__pycache__")), [])


class TestDrugDiscoveryPureEngineAndE2E(unittest.TestCase):
    def test_0_common_graph_alerts_and_canonicalization(self) -> None:
        sys.path.insert(0, str(SCRIPTS_DIR))
        try:
            import _common  # noqa: PLC0415

            # Canonical equivalence & salt stripping
            s1 = _common.standardize_smiles("OCC")
            s2 = _common.standardize_smiles("CCO.Cl.[Na+]")
            self.assertTrue(s1["is_valid"] and s2["is_valid"])
            self.assertEqual(s1["structure_key"], s2["structure_key"])
            self.assertTrue(s2["had_salt"])

            # Neutralization of carboxylate salt
            acet = _common.standardize_smiles("CC(=O)[O-].[Na+]")
            self.assertTrue(acet["is_valid"])
            self.assertTrue(acet["neutralized"])

            # Invalid SMILES detection
            for bad in ("", "C(C", "C1CC", "XYZ123", "C1=CC1(C)(C)(C)C"):
                self.assertFalse(_common.standardize_smiles(bad)["is_valid"], f"expected invalid: {bad}")

            # 15+ PAINS/Brenk alerts coverage check
            alert_samples = {
                "nitroaromatic": "O=[N+]([O-])c1ccc(O)cc1",
                "michael_acceptor": "C=CC(=O)N1CCNCC1",
                "aldehyde": "O=Cc1ccccc1",
                "acyl_halide": "O=C(Cl)c1ccccc1",
                "sulfonyl_halide": "O=S(=O)(Cl)c1ccccc1",
                "epoxide": "CC1CO1",
                "aziridine": "CC1CN1",
                "catechol": "Oc1ccccc1O",
                "quinone": "O=C1C=CC(=O)C=C1",
                "rhodanine": "O=C1NC(=S)SC1=Cc2ccccc2",
                "hydrazine": "NNc1ccccc1",
                "azide": "[N-]=[N+]=Nc1ccccc1",
                "thiourea": "NC(=S)Nc1ccccc1",
                "peroxide": "CC(C)(C)OOC(C)(C)C",
                "alkyl_halide_reactive": "BrCCc1ccccc1",
            }
            for expected_alert, smi in alert_samples.items():
                feat = _common.featurize_smiles(smi)
                self.assertTrue(feat["is_valid"], f"failed to parse {expected_alert}: {smi}")
                self.assertGreater(feat["n_alerts"], 0, f"expected alert for {expected_alert} ({smi})")
        finally:
            if str(SCRIPTS_DIR) in sys.path:
                sys.path.remove(str(SCRIPTS_DIR))

    def test_1_chem_standardize_tsv_and_smi(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            in_tsv = tdir / "raw_compounds.tsv"
            rows = [
                {"compound_id": "CMPD_ASP1", "smiles": "CC(=O)Oc1ccccc1C(=O)O"},
                {"compound_id": "CMPD_ASP_SALT", "smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]"},
                {"compound_id": "CMPD_ETH1", "smiles": "OCC"},
                {"compound_id": "CMPD_ETH2", "smiles": "CCO.Cl"},
                {"compound_id": "CMPD_CHIRAL", "smiles": "C[C@@H](O)c1ccccc1"},
                {"compound_id": "CMPD_BAD", "smiles": "C(C1CC"},
            ]
            with open(in_tsv, "wt", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["compound_id", "smiles"], delimiter="\t")
                w.writeheader()
                w.writerows(rows)

            res = _run_json(
                "chem_standardize.py",
                ["--compounds", str(in_tsv), "--output_dir", str(tdir / "std_out")],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "standardize-result.json")
            self.assertEqual(res["n_input"], 6)
            self.assertEqual(res["n_valid"], 5)
            self.assertEqual(res["n_invalid"], 1)
            self.assertEqual(res["n_salts_stripped"], 2)
            self.assertEqual(res["n_duplicates"], 2)
            self.assertEqual(res["n_output"], 3)
            self.assertFalse(res["keep_duplicates"])

            # Verify keep_duplicates=True retains all 5 valid compounds
            res_keep = _run_json(
                "chem_standardize.py",
                ["--compounds", str(in_tsv), "--keep_duplicates", "--output_dir", str(tdir / "std_keep")],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_keep, "standardize-result.json")
            self.assertTrue(res_keep["keep_duplicates"])
            self.assertEqual(res_keep["n_output"], 5)

            # Verify .smi input support
            smi_path = tdir / "lib.smi"
            smi_path.write_text("CC(=O)Nc1ccc(O)cc1 ACETAMINOPHEN\nCn1cnc2c1c(=O)n(C)c(=O)n2C CAFFEINE\n", encoding="utf-8")
            res_smi = _run_json(
                "chem_standardize.py",
                ["--compounds", str(smi_path), "--output_dir", str(tdir / "std_smi")],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_smi, "standardize-result.json")
            self.assertEqual(res_smi["n_output"], 2)

    def test_2_chem_medchem_triage_tiers_alerts_and_butina(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            in_tsv = tdir / "library.tsv"
            rows = [
                # Lead-like clean analogs (same scaffold -> cluster together)
                {"compound_id": "L1", "smiles": "CC(=O)Nc1ccc(OCCC)cc1"},
                {"compound_id": "L2", "smiles": "CC(=O)Nc1ccc(OCCCC)cc1"},
                # Drug-like clean (MW > 350 or higher cLogP)
                {"compound_id": "D1", "smiles": "Cc1ccc(C(=O)Nc2ccc(C)c(Nc3nccc(-c4cccnc4)n3)c2)cc1"},
                # Alert flagged (PAINS catechol / rhodanine & Brenk nitroaromatic)
                {"compound_id": "A1_RHOD", "smiles": "O=C1NC(=S)SC1=Cc2ccccc2"},
                {"compound_id": "A2_NITRO", "smiles": "O=[N+]([O-])c1ccc(Cl)cc1"},
                # Beyond Ro5 (MW > 500, cLogP > 5, free of PAINS/Brenk alerts)
                {"compound_id": "B1_BRO5", "smiles": "Cc1ccc(C(=O)Nc2ccc(C(C)(C)c3ccc(NC(=O)c4ccc(C(F)(F)F)cc4)cc3)cc2)cc1"},
            ]
            with open(in_tsv, "wt", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["compound_id", "smiles"], delimiter="\t")
                w.writeheader()
                w.writerows(rows)

            res = _run_json(
                "chem_medchem_triage.py",
                ["--compounds", str(in_tsv), "--cluster_cutoff", "0.45", "--output_dir", str(tdir / "triage_out")],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "medchem-triage-result.json")
            self.assertEqual(res["n_valid"], 6)
            self.assertEqual(res["n_alert_flagged"], 2)
            self.assertGreaterEqual(res["tier_counts"]["lead_like"], 1)
            self.assertEqual(res["tier_counts"]["alert_flagged"], 2)
            self.assertEqual(res["tier_counts"]["beyond_ro5"], 1)
            self.assertGreaterEqual(res["n_scaffolds"], 3)

            triage_rows = {r["compound_id"]: r for r in _read_tsv(res["outputs"]["medchem_triage"])}
            self.assertEqual(triage_rows["A1_RHOD"]["triage_tier"], "alert_flagged")
            self.assertEqual(triage_rows["B1_BRO5"]["triage_tier"], "beyond_ro5")
            self.assertEqual(triage_rows["L1"]["cluster_id"], triage_rows["L2"]["cluster_id"])

    def test_3_chem_sar_cliff_le_lle_and_sali(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            in_tsv = tdir / "sar.tsv"
            # Conserved phenyl-amide series with a 1000-fold activity cliff (10 nM vs 10000 nM)
            rows = [
                {"compound_id": "SAR_1", "smiles": "Cc1ccc(C(=O)Nc2ccccc2)cc1", "ic50_nm": "10.0"},
                {"compound_id": "SAR_2", "smiles": "Cc1ccc(C(=O)Nc2ccccc2C)cc1", "ic50_nm": "15.0"},
                {"compound_id": "SAR_3", "smiles": "Cc1ccc(C(=O)Nc2ccc(F)cc2)cc1", "ic50_nm": "12000.0"},
                {"compound_id": "SAR_4", "smiles": "COc1ccc(C(=O)Nc2ccc(Cl)cc2)cc1", "ic50_nm": "500.0"},
            ]
            with open(in_tsv, "wt", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["compound_id", "smiles", "ic50_nm"], delimiter="\t")
                w.writeheader()
                w.writerows(rows)

            res = _run_json(
                "chem_sar_cliff.py",
                [
                    "--compounds", str(in_tsv),
                    "--potency_col", "ic50_nm",
                    "--potency_scale", "nM",
                    "--similarity_threshold", "0.50",
                    "--output_dir", str(tdir / "sar_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "sar-cliff-result.json")
            self.assertEqual(res["n_compounds"], 4)
            self.assertEqual(res["n_pairs"], 6)
            self.assertGreaterEqual(res["n_activity_cliffs"], 1)
            self.assertGreater(res["max_sali"], 3.0)

            sar_rows = {r["compound_id"]: r for r in _read_tsv(res["outputs"]["sar_metrics"])}
            self.assertAlmostEqual(float(sar_rows["SAR_1"]["pic50"]), 8.0, places=3)
            self.assertGreater(float(sar_rows["SAR_1"]["le"]), 0.4)

    def test_4_chem_admet_profile_boiled_egg_esol_herg(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            in_tsv = tdir / "admet.tsv"
            rows = [
                # CNS drug-like BBB+ (Diazepam-like / Caffeine-like)
                {"compound_id": "CNS_CAFF", "smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C"},
                # Lipophilic basic amine -> hERG risk + Pfizer 3/75 risk
                {"compound_id": "HERG_RISK", "smiles": "CN(C)CCCC(c1ccc(Cl)cc1)c2ccc(Cl)cc2"},
                # High TPSA poorly absorbed polar peptide-like analog
                {"compound_id": "POLAR_LOW", "smiles": "NC(=O)C(NC(=O)C(N)CC(=O)O)CC(=O)O"},
            ]
            with open(in_tsv, "wt", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["compound_id", "smiles"], delimiter="\t")
                w.writeheader()
                w.writerows(rows)

            res = _run_json(
                "chem_admet_profile.py",
                ["--compounds", str(in_tsv), "--output_dir", str(tdir / "admet_out")],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "admet-profile-result.json")
            self.assertEqual(res["n_valid"], 3)
            self.assertGreaterEqual(res["n_herg_flagged"], 1)
            self.assertGreaterEqual(res["n_pfizer_3_75_flagged"], 1)
            self.assertGreaterEqual(res["boiled_egg_counts"]["poorly_absorbed"], 1)

            admet_rows = {r["compound_id"]: r for r in _read_tsv(res["outputs"]["admet_profile"])}
            self.assertEqual(admet_rows["HERG_RISK"]["herg_risk"], "True")
            self.assertEqual(admet_rows["HERG_RISK"]["pfizer_3_75_risk"], "True")
            self.assertEqual(admet_rows["POLAR_LOW"]["boiled_egg_class"], "Poorly Absorbed")

    def test_5_chem_docking_eval_pains_sile_consensus_and_fail(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            scores_tsv = tdir / "vina_scores.tsv"
            ligs_tsv = tdir / "ligands.tsv"
            lig_rows = [
                # PAINS rhodanine decoy with artificially strong Vina score
                {"compound_id": "DECOY_RHOD", "smiles": "O=C1NC(=S)SC1=Cc2ccc(Cl)cc2"},
                # Clean compact hit with high LE_dock and SILE
                {"compound_id": "HIT_COMPACT", "smiles": "Cc1ccc(C(=O)Nc2ncccn2)cc1"},
                # Close analog of HIT_COMPACT
                {"compound_id": "HIT_ANALOG", "smiles": "CCc1ccc(C(=O)Nc2ncccn2)cc1"},
                # Distinct clean scaffold hit
                {"compound_id": "HIT_INDOL", "smiles": "c1ccc2c(c1)[nH]cc2CCN3CCOCC3"},
            ]
            score_rows = [
                {"compound_id": "DECOY_RHOD", "vina_kcal": "-10.8", "cnn_score": "0.42"},
                {"compound_id": "HIT_COMPACT", "vina_kcal": "-9.1", "cnn_score": "0.89"},
                {"compound_id": "HIT_COMPACT", "vina_kcal": "-8.3", "cnn_score": "0.75"},  # weaker pose
                {"compound_id": "HIT_ANALOG", "vina_kcal": "-8.9", "cnn_score": "0.81"},
                {"compound_id": "HIT_INDOL", "vina_kcal": "-8.6", "cnn_score": "0.86"},
            ]
            with open(ligs_tsv, "wt", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["compound_id", "smiles"], delimiter="\t")
                w.writeheader()
                w.writerows(lig_rows)
            with open(scores_tsv, "wt", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["compound_id", "vina_kcal", "cnn_score"], delimiter="\t")
                w.writeheader()
                w.writerows(score_rows)

            res = _run_json(
                "chem_docking_eval.py",
                [
                    "--docking_scores", str(scores_tsv),
                    "--score_col", "vina_kcal",
                    "--aux_score_col", "cnn_score",
                    "--compounds", str(ligs_tsv),
                    "--top_n", "3",
                    "--output_dir", str(tdir / "dock_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "docking-eval-result.json")
            self.assertEqual(res["n_input_rows"], 5)
            self.assertEqual(res["n_unique_compounds"], 4)
            self.assertEqual(res["n_pains_filtered"], 1)
            self.assertEqual(res["n_selected_hits"], 3)
            self.assertFalse(res["higher_is_better"])
            # Top clean hit should outrank DECOY_RHOD despite DECOY_RHOD having -10.8 kcal/mol
            self.assertNotEqual(res["top_hits"][0]["compound_id"], "DECOY_RHOD")

            # Verify clean JSON failure on missing file
            err = _run_json(
                "chem_standardize.py",
                ["--compounds", str(tdir / "does_not_exist.tsv"), "--output_dir", str(tdir / "err_out")],
                cwd=tdir,
                expect_fail=True,
            )
            self.assertIn("error", err)


if __name__ == "__main__":
    unittest.main()
