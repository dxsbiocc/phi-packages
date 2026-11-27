#!/usr/bin/env python3
"""Repo-level tests for the metabolomics plugin (outside the plugin directory, so
nothing here ships in the plugin archive — package contract §2 takes git-tracked
files at the package root only).

Run:  python3 -m pytest tests/test_metabolomics_plugin.py
  or: python3 -m unittest tests.test_metabolomics_plugin

Structural tests need pyyaml + jsonschema. Functional tests additionally need the
plugin environment's packages (pandas/numpy/scipy/matplotlib) and are skipped with
a reason when those are absent. Fixtures are synthetic and generated deterministically
into a temp dir; no real data or network is used.

Checks mirror docs/contracts/{plugin,skill,environment,agent,artifact}.md and the
metagenomics tests/test_metagenomics_plugin.py precedent.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
PLUGIN = ROOT / "resources/plugins/metabolomics"
ARTIFACT_SCHEMA = json.loads((ROOT / "docs/contracts/artifact.schema.json").read_text())

try:
    import jsonschema
except ImportError:  # pragma: no cover
    jsonschema = None

NAME_RE = re.compile(r"^[a-z0-9][a-z0-9-]{0,63}$")
TOOL_RE = re.compile(r"^[a-z][a-z0-9_]{0,31}$")
ARG_RE = re.compile(r"^[a-z][a-z0-9_]*$")
PATH_FORMATS = {"input-path", "project-path"}
TYPES = {"string", "number", "integer", "boolean", "array"}
PLATFORMS = {"darwin-arm64": "osx-arm64", "darwin-x64": "osx-64", "linux-x64": "linux-64"}


def frontmatter(path: Path):
    text = path.read_text(encoding="utf-8")
    assert text.startswith("---\n"), f"{path} needs frontmatter"
    head, body = text.split("\n---\n", 1)
    return yaml.safe_load(head[4:]), body


def markdown_links(text: str):
    for link in re.findall(r"\[[^\]]+\]\(([^)\s]+)\)", text):
        if not link.startswith(("http://", "https://", "#")):
            yield link.split("#")[0]


def script_tools():
    for sd in sorted((PLUGIN / "skills").iterdir()):
        fm, _ = frontmatter(sd / "SKILL.md")
        for tool in fm.get("phi", {}).get("scripts", []):
            yield sd, tool


@unittest.skipIf(jsonschema is None, "jsonschema not installed")
class PluginStructureTests(unittest.TestCase):
    def test_manifest_fields(self):
        m = yaml.safe_load((PLUGIN / "phi-package.yaml").read_text())
        self.assertEqual(m["schemaVersion"], 1)
        self.assertEqual(m["id"], PLUGIN.name)
        self.assertEqual(m["type"], "plugin")
        self.assertRegex(m["version"], r"^\d+\.\d+\.\d+(-[0-9A-Za-z.-]+)?$")
        self.assertTrue(1 <= len(m["title"]) <= 80)
        self.assertTrue(1 <= len(m["summary"].strip()) <= 300)
        self.assertRegex(m["toolPrefix"], r"^[a-z][a-z0-9]{1,11}$")

    def test_tool_prefix_unique(self):
        m = yaml.safe_load((PLUGIN / "phi-package.yaml").read_text())
        prefix = m["toolPrefix"]
        for other in (ROOT / "resources/plugins").iterdir():
            if other.name == PLUGIN.name or not (other / "phi-package.yaml").exists():
                continue
            om = yaml.safe_load((other / "phi-package.yaml").read_text())
            self.assertNotEqual(om.get("toolPrefix"), prefix, other.name)
        for sk in (ROOT / "resources/skills").iterdir():
            if (sk / "SKILL.md").exists():
                sfm, _ = frontmatter(sk / "SKILL.md")
                self.assertNotEqual(sfm.get("phi", {}).get("toolPrefix"), prefix, sk.name)

    def test_no_reserved_dirs(self):
        for reserved in ("mcp", "wrappers", "orchestrator"):
            self.assertFalse((PLUGIN / reserved).exists(), reserved)

    def test_components_complete(self):
        m = yaml.safe_load((PLUGIN / "phi-package.yaml").read_text())
        agents = set(m["components"]["agents"])
        skills = set(m["components"]["skills"])
        for rel in agents | skills:
            self.assertTrue((PLUGIN / rel).exists(), rel)
        self.assertEqual({str(p.relative_to(PLUGIN)) for p in (PLUGIN / "agents").glob("*.md")}, agents)
        self.assertEqual({str(p.relative_to(PLUGIN)) for p in (PLUGIN / "skills").iterdir() if p.is_dir()}, skills)

    def test_environments_and_locks(self):
        m = yaml.safe_load((PLUGIN / "phi-package.yaml").read_text())
        for name, e in m.get("environments", {}).items():
            spec_path = PLUGIN / e["spec"]
            spec = yaml.safe_load(spec_path.read_text())
            self.assertEqual(spec["name"], name)
            self.assertRegex(spec["name"], r"^[a-z][a-z0-9-]{0,62}$")
            self.assertTrue(set(spec) <= {"name", "channels", "dependencies", "description", "host", "sourcePackages"})
            self.assertTrue(spec["channels"])
            self.assertGreaterEqual(len(spec["dependencies"]), 1)
            for plat, sub in PLATFORMS.items():
                lock = spec_path.parent / "locks" / f"{plat}.txt"
                self.assertTrue(lock.exists(), f"{name}/{plat}")
                lines = [l.strip() for l in lock.read_text().splitlines()
                         if l.strip() and not l.strip().startswith("#")]
                self.assertEqual(lines[0], "@EXPLICIT", plat)
                for l in lines[1:]:
                    self.assertRegex(l, r"^https://\S+#[0-9a-f]{32}$", plat)
                    self.assertTrue(f"/{sub}/" in l or "/noarch/" in l, f"{plat}: {l[:80]}")

    def test_agent_contract(self):
        names = set()
        for other in (ROOT / "resources/plugins").iterdir():
            for af in (other / "agents").glob("*.md"):
                names.add(frontmatter(af)[0]["name"])
        for af in (PLUGIN / "agents").glob("*.md"):
            a, _ = frontmatter(af)
            self.assertRegex(a["name"], r"^[A-Z][A-Za-z0-9]*$")
            self.assertFalse(a["name"].endswith("Agent"))
            self.assertEqual(a["name"], af.stem)
            self.assertTrue(1 <= len(a["description"].strip()) <= 1024)
            self.assertTrue(a["tools"])
            for reserved in ("spawns", "outputSchema"):
                self.assertNotIn(reserved, a)
            self.assertNotEqual(a.get("visibility", "entry"), "internal")
            self.assertIn(a.get("delegationMode", "optional"), ("required-first", "preferred", "optional"))
            if "fallback" in a:
                fb = a["fallback"]
                self.assertIsInstance(fb["afterFailures"], int)
                self.assertGreater(fb["afterFailures"], 0)
            for s in a.get("skills", []):
                self.assertTrue((PLUGIN / "skills" / s / "SKILL.md").exists(), s)
            names.discard(a["name"])
        for af in (PLUGIN / "agents").glob("*.md"):
            self.assertNotIn(frontmatter(af)[0]["name"], names)

    def test_skill_frontmatter_and_tools(self):
        final_names = set()
        for sd, tool in script_tools():
            fm, _ = frontmatter(sd / "SKILL.md")
            self.assertEqual(fm["name"], sd.name)
            self.assertRegex(fm["name"], NAME_RE)
            self.assertTrue(1 <= len(fm["description"].strip()) <= 1024)
            phi = fm["phi"]
            self.assertNotIn("toolPrefix", phi)
            self.assertEqual(phi.get("environment"), "plugin:metabolomics")
            self.assertTrue(set(phi.get("attachTo", ["main"])) <= {"main", "Metabolomics"})
            tn = tool["name"]
            self.assertRegex(tn, TOOL_RE)
            final = f"mb_{tn}"
            self.assertNotIn(final, final_names)
            final_names.add(final)
            self.assertTrue(1 <= len(tool["description"].strip()) <= 1024)
            self.assertEqual(tool["run"][0], "python")
            script = (sd / tool["run"][1]).resolve()
            self.assertTrue(script.is_file())
            self.assertTrue(str(script).startswith(str(sd.resolve())))
            args = tool["args"]
            self.assertEqual(args["type"], "object")
            self.assertIs(args["additionalProperties"], False)
            for pn, ps in args["properties"].items():
                self.assertRegex(pn, ARG_RE)
                self.assertIn(ps["type"], TYPES)
                if ps["type"] == "array":
                    self.assertIn(ps["items"]["type"], TYPES - {"array"})
                if "format" in ps:
                    self.assertIn(ps["format"], PATH_FORMATS)
                    self.assertEqual(ps["type"], "string")
            for r in args.get("required", []):
                self.assertIn(r, args["properties"])
            self.assertIn(tool["approval"], ("read", "write"))
            if tool["approval"] == "read":
                self.assertTrue(all(ps.get("format") != "project-path"
                                    for ps in args["properties"].values()), final)
            self.assertTrue(1 <= tool.get("timeoutSeconds", 600) <= 86400)
            if "output" in tool:
                sp = (sd / tool["output"]).resolve()
                self.assertTrue(sp.is_file())
                jsonschema.Draft7Validator.check_schema(json.loads(sp.read_text()))
        self.assertEqual(len(final_names), 10)

    def test_markdown_links_resolve(self):
        for sd in sorted((PLUGIN / "skills").iterdir()):
            for md in [sd / "SKILL.md"] + list((sd / "references").rglob("*.md")):
                for link in markdown_links(md.read_text()):
                    self.assertTrue((md.parent / link).resolve().exists(),
                                    f"{md.relative_to(PLUGIN)}: {link}")


def _deps_available():
    try:
        import pandas, numpy, scipy, matplotlib  # noqa: F401
        return True
    except ImportError:
        return False


def make_fixtures(root: Path):
    """Synthetic LC-MS inputs, deterministic (seed 7)."""
    inp = root / "in"
    inp.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)
    n_feat = 20
    qc = [f"QC_{i+1}" for i in range(4)]
    case = [f"C_{i+1}" for i in range(6)]
    control = [f"N_{i+1}" for i in range(6)]
    bio = case + control
    feats = [f"F{i+1}" for i in range(n_feat)]
    base = rng.uniform(1e4, 1e6, n_feat)
    X = pd.DataFrame({c: base * rng.normal(1, 0.05, n_feat) for c in qc}, index=feats)
    for c in bio:
        col = base * rng.normal(1, 0.10, n_feat)
        if c in case:
            col[:4] *= 4.0  # F1..F4 differential
        X[c] = col
    # F5: deterministic high QC RSD (~54%) -> removed by RSD
    b5 = float(X.loc["F5", bio].mean())
    X.loc["F5", qc] = [0.5 * b5, 1.0 * b5, 1.8 * b5, 1.0 * b5]
    # F6: detected in only 2/12 biological samples -> removed by detection
    X.loc["F6", bio[2:]] = np.nan
    # F7: one missing value; F8: one zero (both imputable)
    X.loc["F7", "C_2"] = np.nan
    X.loc["F8", "N_3"] = 0.0
    meta = pd.DataFrame({"mz": 100 + np.arange(n_feat) * 10.5,
                         "rt": 1.0 + np.arange(n_feat) * 0.1}, index=feats)
    out = pd.concat([meta, X], axis=1)
    out.index.name = "feature"
    out.to_csv(inp / "features.tsv", sep="\t")
    md = pd.DataFrame({"group": ["case"] * 6 + ["control"] * 6}, index=bio)
    md.index.name = "sample"
    md.to_csv(inp / "metadata.tsv", sep="\t")
    # MetaboAnalyst-style samples-as-rows table (sample ID, label, then features)
    ma = pd.DataFrame({"sample": ["S1", "S2"], "label": ["case", "control"],
                       "F1": [100.0, 25.0], "F2": [200.0, 50.0]})
    ma.to_csv(inp / "metaboanalyst.csv", index=False)
    # annotation fixture: neutral M=300 cluster + unrelated feature (RT in sec)
    ann = pd.DataFrame({
        "mz": [301.007276, 322.989218, 302.010631, 282.996711, 500.0],
        "rt": [100.0, 100.05, 100.02, 100.03, 300.0],
        "S1": [1000] * 5, "S2": [1100] * 5,
    }, index=["A_MH", "B_MNa", "C_13C", "D_H2O", "E_unrel"])
    ann.index.name = "feature"
    ann.to_csv(inp / "ann-features.tsv", sep="\t")
    # MSP library: NIST + MoNA dialects
    (inp / "lib.msp").write_text(
        "Name: CompoundA\nPrecursorMZ: 301.0073\nNum Peaks: 3\n"
        "100.0 999\n150.05 500\n200.1 300\n\n"
        "NAME: CompoundB\nPRECURSORMZ: 301.0073\nNum Peaks: 3\n"
        "110.0 999\n160.05 500\n210.1 300\n\n"
        "Name: CompoundC\nPrecursorMZ: 400.0\nNum Peaks: 2\n"
        "100.0 999\n150.05 500\n")
    (inp / "query.mgf").write_text(
        "BEGIN IONS\nTITLE=A_MH\nPEPMASS=301.0073\nCHARGE=1+\n"
        "100.0 999\n150.05 500\n200.1 300\nEND IONS\n"
        "BEGIN IONS\nTITLE=nohit\nPEPMASS=999.0\n50.0 100\nEND IONS\n")
    # pathway fixtures: 5 pathways x 4 compounds
    comps = [f"C{i}" for i in range(1, 21)]
    mapping = pd.DataFrame([(f"P{i//4+1}", c) for i, c in enumerate(comps)])
    mapping.to_csv(inp / "mapping.tsv", sep="\t", header=False, index=False)
    (inp / "ids.txt").write_text("C1\nC2\nC3\nC4\n")
    db = pd.DataFrame({"compound_id": comps,
                       "mass": [100.0 + 10 * i for i in range(4)]
                       + [200.0 + i for i in range(16)]})
    db.to_csv(inp / "compound-db.tsv", sep="\t", index=False)
    sig = [100.0 + 10 * i + 1.007276 for i in range(4)]
    bgm = sig + [200.0 + i + 1.007276 for i in range(16)]
    (inp / "sig-mz.txt").write_text("\n".join(f"{m:.6f}" for m in sig) + "\n")
    (inp / "bg-mz.txt").write_text("\n".join(f"{m:.6f}" for m in bgm) + "\n")
    return inp


@unittest.skipIf(jsonschema is None, "jsonschema not installed")
@unittest.skipIf(not _deps_available(), "pandas/numpy/scipy/matplotlib not installed")
class ScriptFunctionalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory(prefix="mb-test-")
        cls.root = Path(cls._tmp.name)
        cls.inp = make_fixtures(cls.root)
        cls.out = cls.root / "out"
        cls.out.mkdir()
        cls.tools = {}
        for sd, tool in script_tools():
            cls.tools[tool["name"]] = (sd, tool)
        cls.artifact_validator = jsonschema.Draft202012Validator(ARTIFACT_SCHEMA)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def script(self, name):
        sd, tool = self.tools[name]
        return (sd / tool["run"][1]).resolve()

    def run_tool(self, name, argv, expect_ok=True):
        sd, tool = self.tools[name]
        r = subprocess.run([sys.executable, str(self.script(name))] + argv,
                           capture_output=True, text=True, cwd=str(self.root))
        if not expect_ok:
            self.assertNotEqual(r.returncode, 0, f"mb_{name} should fail: {r.stdout[:200]}")
            return r
        self.assertEqual(r.returncode, 0, f"mb_{name} failed: {r.stderr[-400:]}")
        obj = json.loads(r.stdout)
        if "output" in tool:
            schema = json.loads((sd / tool["output"]).read_text())
            jsonschema.Draft7Validator(schema).validate(obj)
        self.check_artifacts(obj, name, sd.name)
        return obj

    def check_artifacts(self, obj, name, skill):
        for path in obj.get("artifacts", []):
            p = Path(path)
            self.assertTrue(p.is_file(), f"artifact missing: {path}")
            desc = p.with_name(p.name + ".phi-artifact.json")
            self.assertTrue(desc.is_file(), f"descriptor missing: {desc}")
            d = json.loads(desc.read_text())
            errs = list(self.artifact_validator.iter_errors(d))
            self.assertFalse(errs, f"{desc}: {errs[0].message if errs else ''}")
            self.assertEqual(d["file"], p.name)
            self.assertEqual(d["provenance"]["tool"], f"mb_{name}")
            self.assertEqual(d["provenance"]["skill"], skill)

    def o(self, *p):
        return str(self.out.joinpath(*p))

    def i(self, name):
        return str(self.inp / name)

    def filtered(self):
        path = self.o("qc", "filtered-matrix.tsv")
        if not Path(path).exists():
            self.run_tool("qc", ["--input", self.i("features.tsv"),
                                 "--output_dir", self.o("qc")])
        return path

    def processed(self):
        path = self.o("processed.tsv")
        if not Path(path).exists():
            self.run_tool("preprocess", [
                "--input", self.filtered(), "--output", path,
                "--impute", "half_min", "--normalize", "pqn", "--transform", "log2"])
        return path

    def processed_bio(self):
        path = self.o("processed-bio.tsv")
        if not Path(path).exists():
            m = pd.read_csv(self.processed(), sep="\t", index_col=0)
            bio = [c for c in m.columns if c.startswith(("C_", "N_"))]
            m[bio].to_csv(path, sep="\t")
        return path

    def test_flag_parity(self):
        for name, (sd, tool) in self.tools.items():
            r = subprocess.run([sys.executable, str(self.script(name)), "--help"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, f"mb_{name} --help: {r.stderr[-200:]}")
            flags = {f.lstrip("-") for f in re.findall(r"(?<![\w-])(--[a-z][a-z0-9_-]*)", r.stdout)} - {"help"}
            missing = set(tool["args"]["properties"]) - flags
            self.assertFalse(missing, f"mb_{name}: declared args not flags: {missing}")

    def test_inspect(self):
        obj = self.run_tool("inspect", ["--input", self.i("features.tsv")])
        f = obj["files"][0]
        self.assertEqual(f["kind"], "feature_table")
        self.assertEqual(f["format"], "xcms")
        self.assertEqual(f["n_features"], 20)
        self.assertEqual(f["n_samples"], 16)
        self.assertEqual(f["rt_unit"], "min")
        self.assertEqual(len(f["qc_samples"]), 4)
        self.assertIn("qc", f["ready_steps"])
        self.assertIn("annotate", f["ready_steps"])

    def test_inspect_other_kinds(self):
        obj = self.run_tool("inspect", ["--input", str(self.inp)])
        kinds = {Path(f["path"]).name: f["kind"] for f in obj["files"]}
        self.assertEqual(kinds["metaboanalyst.csv"], "feature_table")
        ma = next(f for f in obj["files"] if f["path"].endswith("metaboanalyst.csv"))
        self.assertEqual(ma["format"], "metaboanalyst")
        self.assertEqual(ma["orientation"], "samples-as-rows")
        self.assertEqual(kinds["metadata.tsv"], "metadata")
        self.assertEqual(kinds["lib.msp"], "msp_library")
        self.assertEqual(kinds["query.mgf"], "mgf_spectra")
        self.assertEqual(kinds["mapping.tsv"], "pathway_mapping")
        self.assertEqual(kinds["compound-db.tsv"], "compound_db")
        self.assertEqual(kinds["ids.txt"], "id_list")
        self.assertEqual(kinds["sig-mz.txt"], "mass_list")

    def test_qc_filtering(self):
        obj = self.run_tool("qc", ["--input", self.i("features.tsv"),
                                   "--output_dir", self.o("qc")])
        self.assertEqual(obj["n_features_before"], 20)
        self.assertEqual(obj["removed_by_rsd"], 1)      # F5
        self.assertEqual(obj["removed_by_detection"], 1)  # F6
        self.assertEqual(obj["n_features_after"], 18)
        report = pd.read_csv(self.o("qc", "qc-report.tsv"), sep="\t", index_col=0)
        self.assertFalse(report.loc["F5", "kept"])
        self.assertFalse(report.loc["F6", "kept"])
        self.assertTrue(report.loc["F1", "kept"])

    def test_qc_requires_qc_samples(self):
        df = pd.read_csv(self.i("features.tsv"), sep="\t", index_col=0)
        bio = [c for c in df.columns if not c.startswith("QC")]
        noqc = self.o("noqc.tsv")
        df[bio].to_csv(noqc, sep="\t")
        r = self.run_tool("qc", ["--input", noqc, "--output_dir", self.o("noqc-out")],
                          expect_ok=False)
        self.assertIn("no QC samples", r.stdout)

    def test_qc_determinism(self):
        a = self.o("qc-det-a")
        b = self.o("qc-det-b")
        self.run_tool("qc", ["--input", self.i("features.tsv"), "--output_dir", a])
        self.run_tool("qc", ["--input", self.i("features.tsv"), "--output_dir", b])
        for name in ("filtered-matrix.tsv", "qc-report.tsv"):
            self.assertEqual(Path(a, name).read_bytes(), Path(b, name).read_bytes())

    def test_preprocess_chain(self):
        obj = self.run_tool("preprocess", [
            "--input", self.filtered(), "--output", self.o("processed.tsv"),
            "--impute", "half_min", "--normalize", "pqn", "--transform", "log2"])
        self.assertEqual(obj["n_imputed"], 2)  # F7 NaN + F8 zero (F5/F6 already filtered)
        self.assertEqual(obj["n_missing_remaining"], 0)
        m = pd.read_csv(self.o("processed.tsv"), sep="\t", index_col=0)
        self.assertTrue((m.select_dtypes("number") > -30).all().all())
        factors = obj["factors"]
        self.assertTrue(all(0.5 < v < 2.0 for v in factors.values()))

    def test_preprocess_log_rejects_zeros(self):
        r = self.run_tool("preprocess", [
            "--input", self.i("features.tsv"), "--output", self.o("bad.tsv"),
            "--transform", "log2"], expect_ok=False)
        self.assertIn("error", r.stdout)

    def test_preprocess_tic_and_is(self):
        obj = self.run_tool("preprocess", [
            "--input", self.filtered(), "--output", self.o("tic.tsv"),
            "--impute", "half_min", "--normalize", "tic"])
        m = pd.read_csv(self.o("tic.tsv"), sep="\t", index_col=0)
        sums = m[[c for c in m.columns if c not in ("mz", "rt")]].sum()
        self.assertTrue(np.allclose(sums, 1.0))
        # internal standard: add an IS column, normalize by it, column dropped
        df = pd.read_csv(self.filtered(), sep="\t", index_col=0)
        df["IS"] = df[[c for c in df.columns if c.startswith("C_")]].mean(axis=1)
        is_path = self.o("with-is.tsv")
        df.to_csv(is_path, sep="\t")
        obj = self.run_tool("preprocess", [
            "--input", is_path, "--output", self.o("is-norm.tsv"),
            "--impute", "half_min", "--normalize", "is", "--is_col", "IS"])
        self.assertEqual(obj["dropped_columns"], ["IS"])
        m = pd.read_csv(self.o("is-norm.tsv"), sep="\t", index_col=0)
        self.assertNotIn("IS", m.columns)

    def test_multivariate(self):
        obj = self.run_tool("multivariate", [
            "--input", self.processed_bio(), "--metadata", self.i("metadata.tsv"),
            "--output_dir", self.o("mv"), "--n_perm", "19"])
        self.assertEqual(obj["n_samples"], 12)
        self.assertEqual(obj["cv_scheme"], "loo")
        self.assertEqual(obj["groups"], {"case": 6, "control": 6})
        vip_feats = [v["feature"] for v in obj["vip_top"]]
        self.assertTrue(set(vip_feats[:4]) & {"F1", "F2", "F3", "F4"})
        self.assertIn("q2_p", obj["plsda"]["permutation"])
        self.assertTrue(Path(self.o("mv", "pca-scores.png")).is_file())
        self.assertTrue(Path(self.o("mv", "plsda-scores.png")).is_file())

    def test_multivariate_determinism_and_noperm(self):
        a = self.run_tool("multivariate", [
            "--input", self.processed_bio(), "--metadata", self.i("metadata.tsv"),
            "--output_dir", self.o("mv-a"), "--n_perm", "9", "--seed", "42"])
        b = self.run_tool("multivariate", [
            "--input", self.processed_bio(), "--metadata", self.i("metadata.tsv"),
            "--output_dir", self.o("mv-b"), "--n_perm", "9", "--seed", "42"])
        self.assertEqual(a["plsda"]["permutation"]["q2_p"],
                         b["plsda"]["permutation"]["q2_p"])
        self.assertEqual(a["vip_top"], b["vip_top"])
        c = self.run_tool("multivariate", [
            "--input", self.processed_bio(), "--metadata", self.i("metadata.tsv"),
            "--output_dir", self.o("mv-c"), "--n_perm", "0", "--scaling", "pareto"])
        self.assertEqual(c["plsda"]["permutation"], {"n_perm": 0})

    def test_univariate_welch_logged(self):
        obj = self.run_tool("univariate", [
            "--input", self.processed_bio(), "--metadata", self.i("metadata.tsv"),
            "--case", "case", "--control", "control", "--test", "welch",
            "--logged", "--output_dir", self.o("uv-welch")])
        self.assertGreaterEqual(obj["n_significant"], 4)
        top = {t["feature"]: t for t in obj["top"]}
        self.assertIn("F1", top)
        self.assertTrue(1.5 < top["F1"]["log2_fold_change"] < 2.5)
        self.assertTrue(Path(self.o("uv-welch", "volcano.png")).is_file())

    def test_univariate_wilcoxon_and_kruskal(self):
        obj = self.run_tool("univariate", [
            "--input", self.processed_bio(), "--metadata", self.i("metadata.tsv"),
            "--case", "case", "--control", "control", "--test", "wilcoxon",
            "--logged", "--output_dir", self.o("uv-wilcox")])
        self.assertGreaterEqual(obj["n_significant"], 1)
        obj = self.run_tool("univariate", [
            "--input", self.processed_bio(), "--metadata", self.i("metadata.tsv"),
            "--test", "kruskal", "--logged", "--output_dir", self.o("uv-kw")])
        self.assertEqual(obj["test"], "kruskal")
        self.assertFalse(Path(self.o("uv-kw", "volcano.png")).exists())

    def test_univariate_requires_case_control(self):
        r = self.run_tool("univariate", [
            "--input", self.processed_bio(), "--metadata", self.i("metadata.tsv"),
            "--test", "welch", "--output_dir", self.o("uv-bad")], expect_ok=False)
        self.assertIn("case", r.stdout)

    def test_annotate(self):
        obj = self.run_tool("annotate", [
            "--input", self.i("ann-features.tsv"), "--output_dir", self.o("ann")])
        self.assertEqual(obj["rt_unit"], "sec")
        self.assertEqual(obj["rt_tol"], 12.0)
        self.assertEqual(obj["relationships_by_type"],
                         {"adduct_pair": 1, "isotope_13C": 1, "neutral_loss": 1})
        self.assertEqual(obj["n_pseudospectra"], 1)
        grp = pd.read_csv(self.o("ann", "pseudospectra.tsv"), sep="\t")
        self.assertAlmostEqual(grp["neutral_mass"][0], 300.0, places=5)
        self.assertEqual(set(grp["members"][0].split(";")),
                         {"A_MH", "B_MNa", "C_13C", "D_H2O"})
        self.assertIn("level 3", obj["caveat"])

    def test_annotate_neg_mode(self):
        df = pd.read_csv(self.i("ann-features.tsv"), sep="\t", index_col=0)
        neg = pd.DataFrame({
            "mz": [299.0 - 1.007276, 299.0 + 34.969402],
            "rt": [2.0, 2.05], "S1": [1000, 900], "S2": [1100, 950],
        }, index=["M_H", "M_Cl"])
        neg.index.name = "feature"
        p = self.o("neg-features.tsv")
        neg.to_csv(p, sep="\t")
        obj = self.run_tool("annotate", ["--input", p, "--output_dir", self.o("ann-neg"),
                                         "--mode", "neg"])
        self.assertEqual(obj["relationships_by_type"].get("adduct_pair"), 1)
        grp = pd.read_csv(self.o("ann-neg", "pseudospectra.tsv"), sep="\t")
        self.assertAlmostEqual(grp["neutral_mass"][0], 299.0, places=5)

    def test_ms2_match(self):
        obj = self.run_tool("ms2_match", [
            "--query", self.i("query.mgf"), "--library", self.i("lib.msp"),
            "--output_dir", self.o("ms2")])
        self.assertEqual(obj["n_library_spectra"], 3)
        self.assertEqual(obj["n_queries_with_hits"], 1)
        self.assertEqual(obj["n_hits"], 1)
        t = pd.read_csv(self.o("ms2", "ms2-matches.tsv"), sep="\t")
        self.assertEqual(t["library_name"][0], "CompoundA")
        self.assertAlmostEqual(t["cosine"][0], 1.0, places=6)
        self.assertEqual(t["matched_fragments"][0], 3)
        self.assertIn("level 2", obj["caveat"])

    def test_pathway_ids(self):
        obj = self.run_tool("pathway", [
            "--mode", "ids", "--input", self.i("ids.txt"),
            "--mapping", self.i("mapping.tsv"), "--output_dir", self.o("pw-ids")])
        self.assertEqual(obj["background"], {"kind": "compounds", "n": 20})
        top = obj["top"][0]
        self.assertEqual(top["pathway"], "P1")
        self.assertEqual(top["overlap"], 4)
        self.assertAlmostEqual(top["pvalue"], 1.0 / 4845.0, places=8)
        self.assertEqual(obj["n_significant"], 1)

    def test_pathway_masses_and_determinism(self):
        argv = ["--mode", "masses", "--input", self.i("sig-mz.txt"),
                "--mapping", self.i("mapping.tsv"),
                "--compound_db", self.i("compound-db.tsv"),
                "--background", self.i("bg-mz.txt"), "--n_perm", "199"]
        a = self.run_tool("pathway", argv + ["--output_dir", self.o("pw-m-a")])
        b = self.run_tool("pathway", argv + ["--output_dir", self.o("pw-m-b")])
        self.assertEqual(a["top"], b["top"])
        self.assertEqual(a["n_matched_input"], 4)
        self.assertEqual(a["top"][0]["pathway"], "P1")
        self.assertEqual(a["top"][0]["overlap"], 4)
        self.assertLessEqual(a["top"][0]["pvalue"], 0.01)
        self.assertEqual(a["n_significant"], 1)

    def test_pathway_requires_mapping(self):
        r = self.run_tool("pathway", [
            "--mode", "ids", "--input", self.i("ids.txt"),
            "--output_dir", self.o("pw-bad")], expect_ok=False)
        self.assertIn("mapping", r.stdout)

    def test_plot_heatmap(self):
        obj = self.run_tool("plot", [
            "--input", self.processed_bio(), "--output_dir", self.o("hm"),
            "--trans", "log10", "--cluster", "both", "--top", "8"])
        self.assertEqual(obj["n_features_shown"], 8)
        self.assertEqual(obj["n_samples"], 12)
        self.assertTrue(Path(self.o("hm", "feature-heatmap.png")).is_file())

    def test_drift_correct(self):
        # 24 injections (6 QC interleaved with 18 biological samples) with
        # strong multiplicative signal attenuation along injection order (~60% drop)
        rng = np.random.default_rng(2026)
        n_inj = 24
        orders = np.arange(1, n_inj + 1, dtype=float)
        qc_indices = {0, 4, 9, 14, 18, 23}
        sample_names = []
        sample_types = []
        for idx in range(n_inj):
            if idx in qc_indices:
                sample_names.append(f"QC_{idx+1}")
                sample_types.append("QC")
            else:
                grp = "Case" if idx % 2 == 0 else "Ctrl"
                sample_names.append(f"{grp}_{idx+1}")
                sample_types.append(grp)

        sheet_path = Path(self.o("drift_sample_sheet.tsv"))
        sheet_path.parent.mkdir(parents=True, exist_ok=True)
        sheet_lines = ["sample\tinjection_order\tsample_type"]
        for s_name, ord_v, st in zip(sample_names, orders, sample_types):
            sheet_lines.append(f"{s_name}\t{int(ord_v)}\t{st}")
        sheet_path.write_text("\n".join(sheet_lines) + "\n")

        drift_factor = 1.0 - 0.58 * ((orders - 1.0) / (n_inj - 1.0)) + 0.08 * np.sin(orders / 3.5)
        rows = ["feature\tmz\trt\t" + "\t".join(sample_names)]
        for f_i in range(25):
            base = 1e5 * (1.0 + 0.25 * (f_i % 5))
            vals = []
            for j in range(n_inj):
                bio_mult = 1.35 if (sample_types[j] == "Case" and f_i < 6) else 1.0
                noise = rng.normal(1.0, 0.025)
                vals.append(f"{base * bio_mult * drift_factor[j] * noise:.4f}")
            rows.append(f"feat_{f_i:02d}\t{150.0 + f_i * 12.5:.4f}\t{60.0 + f_i * 8.0:.2f}\t" + "\t".join(vals))

        mat_path = Path(self.o("drift_raw_matrix.tsv"))
        mat_path.write_text("\n".join(rows) + "\n")

        out_dir = self.o("drift_out")
        obj = self.run_tool("drift_correct", [
            "--input", str(mat_path),
            "--sample_sheet", str(sheet_path),
            "--output_dir", out_dir,
            "--qc_label", "QC",
            "--span", "0.75",
            "--rsd_cutoff", "15.0",
        ])
        self.assertEqual(obj["tool"], "mb_drift_correct")
        self.assertEqual(obj["n_features"], 25)
        self.assertEqual(obj["n_qc_samples"], 6)
        self.assertLess(obj["median_qc_rsd_after"], 0.5 * obj["median_qc_rsd_before"])
        self.assertGreater(obj["n_features_rescued"], 0)
        corr_mat_path = Path(obj["outputs"]["drift_corrected_matrix"])
        self.assertTrue(corr_mat_path.is_file())
        self.assertTrue(Path(obj["outputs"]["feature_qc_metrics"]).is_file())
        self.assertTrue(Path(obj["outputs"]["drift_qc_plot"]).is_file())
        hdr = corr_mat_path.read_text().splitlines()[0].split("\t")
        self.assertIn("mz", hdr)
        self.assertIn("rt", hdr)
        sc_ctx = obj.get("scientificContext", {})
        self.assertIn(str(corr_mat_path), sc_ctx)
        self.assertEqual(sc_ctx[str(corr_mat_path)]["valueSemantics"], "normalized")
        for art in obj["artifacts"]:
            ap = Path(art)
            self.assertTrue(ap.is_file())
            sidecar = ap.with_name(ap.name + ".phi-artifact.json")
            self.assertTrue(sidecar.is_file())
            meta = json.loads(sidecar.read_text())
            self.assertEqual(meta["contractVersion"], "1.0.0")
            self.assertEqual(meta["provenance"]["tool"], "mb_drift_correct")
            self.assertEqual(meta["provenance"]["skill"], "lcms-workflow")


if __name__ == "__main__":
    unittest.main()