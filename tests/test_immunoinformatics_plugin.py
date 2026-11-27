#!/usr/bin/env python3
"""Contract, layout, AST, schema, and end-to-end runtime tests for immunoinformatics (0.1.0)."""

from __future__ import annotations

import ast
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
PLUGIN_DIR = ROOT / "resources" / "plugins" / "immunoinformatics"
SKILL_DIR = PLUGIN_DIR / "skills" / "immuno-workflow"

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
        testcase.assertEqual(meta["provenance"]["skill"], "immuno-workflow")


class TestImmunoinformaticsContract(unittest.TestCase):
    def test_plugin_layout_locks_schemas_and_ast(self) -> None:
        manifest = (PLUGIN_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("schemaVersion: 1", manifest)
        self.assertIn("id: immunoinformatics", manifest)
        self.assertIn("type: plugin", manifest)
        self.assertIn("version: 0.1.0", manifest)
        self.assertIn("toolPrefix: imm", manifest)
        self.assertIn("- agents/Immunoinformatics.md", manifest)
        self.assertIn("- skills/immuno-workflow", manifest)
        self.assertIn("spec: environments/imm/environment.yml", manifest)

        # Check summary <= 260 characters
        m_sum = re.search(r"^summary:\s*>-\n((?:  .*\n)+)", manifest, re.M)
        self.assertIsNotNone(m_sum)
        summary_text = " ".join(line.strip() for line in m_sum.group(1).splitlines())
        self.assertLessEqual(len(summary_text), 260, f"summary too long: {len(summary_text)}")

        for req in ("icon.svg", "README.md", "CHANGELOG.md", "agents/Immunoinformatics.md"):
            self.assertTrue((PLUGIN_DIR / req).is_file(), f"missing {req}")
        for reserved in ("mcp", "wrappers", "orchestrator"):
            self.assertFalse((PLUGIN_DIR / reserved).exists(), f"reserved dir {reserved} exists")

        agent_txt = (PLUGIN_DIR / "agents/Immunoinformatics.md").read_text(encoding="utf-8")
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
            lock = PLUGIN_DIR / f"environments/imm/locks/{plat}.txt"
            self.assertTrue(lock.is_file(), f"missing lock {plat}")
            ltxt = lock.read_text(encoding="utf-8")
            self.assertIn("@EXPLICIT", ltxt)
            self.assertIn("# spec: resources/plugins/immunoinformatics/environments/imm/environment.yml", ltxt)

        stext = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
        sfm = stext[4 : stext.index("\n---\n", 4)]
        sbody = stext[stext.index("\n---\n", 4) + 5 :]
        self.assertIn("name: immuno-workflow", sfm)
        self.assertIn("  environment: plugin:imm", sfm)
        self.assertIn("    - Immunoinformatics", sfm)
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
                "repertoire_inspect",
                "diversity",
                "clonotype_tracking",
                "hla_summary",
                "neoantigen_prioritize",
            ],
        )
        self.assertEqual(len(outs), 5)
        for rel_out in outs:
            schema_file = SKILL_DIR / rel_out
            self.assertTrue(schema_file.is_file(), f"missing schema {rel_out}")
            json.loads(schema_file.read_text(encoding="utf-8"))

        exp_scripts = {
            "_common.py",
            "imm_repertoire_inspect.py",
            "imm_diversity.py",
            "imm_clonotype_tracking.py",
            "imm_hla_summary.py",
            "imm_neoantigen_prioritize.py",
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


@unittest.skipUnless(HAS_SCIPY_STACK, "numpy/pandas/scipy/sklearn/matplotlib required for live E2E")
class TestImmunoinformaticsRuntimeE2E(unittest.TestCase):
    def test_1_repertoire_inspect_multi_format(self) -> None:
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # 1) AIRR TSV
            airr_path = tdir / "sample_airr.tsv"
            pd.DataFrame(
                [
                    {
                        "sample_id": "S_AIRR",
                        "locus": "TRB",
                        "v_call": "TRBV12-3*01",
                        "j_call": "TRBJ2-7*01",
                        "junction_aa": "CASSLGQAYEQYF",
                        "productive": True,
                        "duplicate_count": 40,
                    },
                    {
                        "sample_id": "S_AIRR",
                        "locus": "TRB",
                        "v_call": "TRBV5-1*01",
                        "j_call": "TRBJ1-1*01",
                        "junction_aa": "CASSL*QAYEQYF",
                        "productive": False,
                        "duplicate_count": 5,
                    },
                    {
                        "sample_id": "S_AIRR",
                        "locus": "TRA",
                        "v_call": "TRAV12-1*01",
                        "j_call": "TRAJ20*01",
                        "junction_aa": "CAVNNDYKLSF",
                        "productive": True,
                        "duplicate_count": 15,
                    },
                ]
            ).to_csv(airr_path, sep="\t", index=False)

            # 2) 10x Cell Ranger VDJ CSV
            cr_path = tdir / "Pt1_filtered_contig_annotations.csv"
            pd.DataFrame(
                [
                    {
                        "barcode": "AAACCTGAGAA-1",
                        "contig_id": "contig_1",
                        "chain": "TRB",
                        "v_gene": "TRBV7-2",
                        "d_gene": "TRBD1",
                        "j_gene": "TRBJ2-1",
                        "c_gene": "TRBC2",
                        "cdr3": "CASSLAPGATNEKLFF",
                        "cdr3_nt": "TGTGCCAGCAGCTTAGCTCCGGGAGCGACTAATGAAAAACTGTTTTTT",
                        "productive": "True",
                        "umis": 12,
                        "reads": 120,
                    },
                    {
                        "barcode": "AAACCTGAGAA-2",
                        "contig_id": "contig_2",
                        "chain": "TRB",
                        "v_gene": "TRBV20-1",
                        "d_gene": "None",
                        "j_gene": "TRBJ1-5",
                        "c_gene": "TRBC1",
                        "cdr3": "CSARDLTSGANVLTF",
                        "cdr3_nt": "TGCAGTGCTAGAGACCTAACCAGCGGGGCAAATGTCCTGACTTTC",
                        "productive": "True",
                        "umis": 8,
                        "reads": 80,
                    },
                ]
            ).to_csv(cr_path, index=False)

            # 3) TRUST4 report TSV (with #count header!)
            trust_path = tdir / "Pt2_report.tsv"
            trust_path.write_text(
                "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\tcid\n"
                "50\t0.5\tTGTGCCAGCAGTTAC\tCASSYSGQGYEQYF\tTRBV6-5*01\tTRBD1\tTRBJ2-7*01\tTRBC2\tcid1\n"
                "30\t0.3\tTGTGCCAGCAGCTTA\tout_of_frame\tTRBV12-3*01\t.\tTRBJ1-1*01\tTRBC1\tcid2\n"
                "20\t0.2\tTGTGCGAGAGGGACT\tCARGTYYYYMDVW\tIGHV3-23*01\tIGHD3-10\tIGHJ6*02\tIGHG1\tcid3\n",
                encoding="utf-8",
            )

            # 4) MiXCR clones TSV
            mixcr_path = tdir / "Pt3_clones.tsv"
            pd.DataFrame(
                [
                    {
                        "cloneId": 0,
                        "cloneCount": 100,
                        "cloneFraction": 0.8,
                        "nSeqCDR3": "TGTGCCAGCAGCTTAGGACAGGGTTACGAGCAGTACTTC",
                        "aaSeqCDR3": "CASSLGQGYEQYF",
                        "allVHitsWithScore": "TRBV12-3*00(150.0),TRBV12-4*00(130.0)",
                        "allJHitsWithScore": "TRBJ2-7*00(90.0)",
                    },
                    {
                        "cloneId": 1,
                        "cloneCount": 25,
                        "cloneFraction": 0.2,
                        "nSeqCDR3": "TGTGCCAGCAGCTTAGGACAGGGTTACGAGCAGTACTTCA",
                        "aaSeqCDR3": "CASSLGQ_GYEQYF",
                        "allVHitsWithScore": "TRBV12-3*00(140.0)",
                        "allJHitsWithScore": "TRBJ2-7*00(85.0)",
                    },
                ]
            ).to_csv(mixcr_path, sep="\t", index=False)

            res = _run_json(
                "imm_repertoire_inspect.py",
                [
                    "--inputs",
                    str(airr_path),
                    str(cr_path),
                    str(trust_path),
                    str(mixcr_path),
                    "--chain_filter",
                    "TCR",
                    "--productive_only",
                    "--output_dir",
                    str(tdir / "rep_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "repertoire-inspect-result.json")
            self.assertEqual(res["n_samples"], 4)
            self.assertEqual(res["detected_formats"]["sample_airr.tsv"], "airr")
            self.assertEqual(res["detected_formats"]["Pt1_filtered_contig_annotations.csv"], "cellranger_vdj")
            self.assertEqual(res["detected_formats"]["Pt2_report.tsv"], "trust4")
            self.assertEqual(res["detected_formats"]["Pt3_clones.tsv"], "mixcr")

            harm = pd.read_csv(res["outputs"]["harmonized_repertoire"], sep="\t")
            self.assertTrue(bool(harm["productive"].all()))
            self.assertTrue((harm["receptor_type"] == "TCR").all())

            qc_df = pd.read_csv(res["outputs"]["repertoire_qc_summary"], sep="\t").set_index("sample_id")
            self.assertGreater(float(qc_df.loc["S_AIRR", "stop_codon_pct"]), 0.0)
            self.assertGreater(float(qc_df.loc["Pt2", "out_of_frame_pct"]), 0.0)

    def test_2_diversity_homeostasis_and_group_stats(self) -> None:
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            rows = []
            # 3 Healthy polyclonal samples (uniform counts) vs 3 Responder oligoclonal samples (1 hyperexpanded clone)
            for idx in range(1, 4):
                sid = f"Healthy_{idx}"
                for k in range(30):
                    rows.append(
                        {
                            "sample_id": sid,
                            "locus": "TRB",
                            "v_call": f"TRBV{(k % 6) + 1}*01",
                            "j_call": f"TRBJ1-{(k % 3) + 1}*01",
                            "cdr3_aa": f"CASSL{chr(65 + (k % 20))}GQAYEQYF",
                            "productive": True,
                            "duplicate_count": 10 + (k % 3),
                        }
                    )
            for idx in range(1, 4):
                sid = f"Resp_{idx}"
                rows.append(
                    {
                        "sample_id": sid,
                        "locus": "TRB",
                        "v_call": "TRBV12-3*01",
                        "j_call": "TRBJ2-7*01",
                        "cdr3_aa": "CASSLGQAYEQYF",
                        "productive": True,
                        "duplicate_count": 800,
                    }
                )
                for k in range(1, 20):
                    rows.append(
                        {
                            "sample_id": sid,
                            "locus": "TRB",
                            "v_call": f"TRBV{(k % 5) + 2}*01",
                            "j_call": "TRBJ1-1*01",
                            "cdr3_aa": f"CASSF{chr(65 + k)}NEKLFF",
                            "productive": True,
                            "duplicate_count": 2,
                        }
                    )
            rep_path = tdir / "rep.tsv"
            pd.DataFrame(rows).to_csv(rep_path, sep="\t", index=False)

            meta_path = tdir / "meta.tsv"
            pd.DataFrame(
                [
                    {"sample_id": f"Healthy_{i}", "response_group": "Healthy"} for i in range(1, 4)
                ]
                + [
                    {"sample_id": f"Resp_{i}", "response_group": "Responder"} for i in range(1, 4)
                ]
            ).to_csv(meta_path, sep="\t", index=False)

            res = _run_json(
                "imm_diversity.py",
                [
                    "--repertoire",
                    str(rep_path),
                    "--metadata",
                    str(meta_path),
                    "--group_col",
                    "response_group",
                    "--locus",
                    "TRB",
                    "--output_dir",
                    str(tdir / "div_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "diversity-result.json")
            self.assertEqual(res["n_samples"], 6)
            div_df = pd.read_csv(res["outputs"]["diversity_metrics"], sep="\t").set_index("sample_id")
            self.assertGreater(float(div_df.loc["Resp_1", "pielou_clonality"]), float(div_df.loc["Healthy_1", "pielou_clonality"]))
            self.assertGreater(float(div_df.loc["Resp_1", "homeostasis_hyperexpanded_pct"]), 80.0)
            tot_homeo = (
                div_df.loc["Healthy_1", "homeostasis_hyperexpanded_pct"]
                + div_df.loc["Healthy_1", "homeostasis_large_pct"]
                + div_df.loc["Healthy_1", "homeostasis_medium_pct"]
                + div_df.loc["Healthy_1", "homeostasis_small_pct"]
            )
            self.assertAlmostEqual(float(tot_homeo), 100.0, delta=0.1)
            self.assertGreater(len(res["group_comparisons"]), 0)
            self.assertIsNotNone(res["outputs"]["diversity_group_comparisons"])

    def test_3_clonotype_tracking_morisita_fisher(self) -> None:
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            rows = [
                # Baseline: 1000 reads total
                {"sample_id": "Pre", "locus": "TRB", "v_call": "TRBV12-3", "cdr3_aa": "CASSLGQAYEQYF", "productive": True, "duplicate_count": 20},
                {"sample_id": "Pre", "locus": "TRB", "v_call": "TRBV7-2", "cdr3_aa": "CASSLAPGATNEKLFF", "productive": True, "duplicate_count": 250},
                {"sample_id": "Pre", "locus": "TRB", "v_call": "TRBV5-1", "cdr3_aa": "CASSFSTCSANYGYTF", "productive": True, "duplicate_count": 100},
                # Followup: 1000 reads total (CASSLGQAYEQYF expands 20 -> 350; CASSLAPGATNEKLFF contracts 250 -> 15)
                {"sample_id": "Post", "locus": "TRB", "v_call": "TRBV12-3", "cdr3_aa": "CASSLGQAYEQYF", "productive": True, "duplicate_count": 350},
                {"sample_id": "Post", "locus": "TRB", "v_call": "TRBV7-2", "cdr3_aa": "CASSLAPGATNEKLFF", "productive": True, "duplicate_count": 15},
                {"sample_id": "Post", "locus": "TRB", "v_call": "TRBV5-1", "cdr3_aa": "CASSFSTCSANYGYTF", "productive": True, "duplicate_count": 100},
            ]
            # Fill background singletons sharing CASSLG motif
            for k in range(30):
                rows.append({"sample_id": "Pre", "locus": "TRB", "v_call": "TRBV20-1", "cdr3_aa": f"CSARD{chr(65 + (k % 20))}NVLTF", "productive": True, "duplicate_count": 21})
                rows.append({"sample_id": "Post", "locus": "TRB", "v_call": "TRBV6-5", "cdr3_aa": f"CASSE{chr(65 + (k % 20))}QETQYF", "productive": True, "duplicate_count": 17})

            rep_path = tdir / "long_rep.tsv"
            pd.DataFrame(rows).to_csv(rep_path, sep="\t", index=False)

            res = _run_json(
                "imm_clonotype_tracking.py",
                [
                    "--repertoire",
                    str(rep_path),
                    "--baseline_sample",
                    "Pre",
                    "--followup_sample",
                    "Post",
                    "--locus",
                    "TRB",
                    "--output_dir",
                    str(tdir / "track_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "clonotype-tracking-result.json")
            self.assertEqual(res["n_samples"], 2)
            self.assertGreaterEqual(res["n_expanded"], 1)
            self.assertGreaterEqual(res["n_contracted"], 1)
            self.assertGreater(res["mean_morisita_horn"], 0.0)
            self.assertLess(res["mean_morisita_horn"], 1.0)

            diff_df = pd.read_csv(res["outputs"]["clonotype_differential_tracking"], sep="\t").set_index("cdr3_aa")
            self.assertEqual(diff_df.loc["CASSLGQAYEQYF", "trend"], "expanded")
            self.assertEqual(diff_df.loc["CASSLAPGATNEKLFF", "trend"], "contracted")

    def test_4_hla_summary_multi_caller_and_hed(self) -> None:
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # 1) OptiType TSV for PtA (heterozygous at A/B/C)
            opt_path = tdir / "PtA_result.tsv"
            opt_path.write_text(
                "\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n"
                "0\tA*01:01\tA*02:01\tB*27:05\tB*44:02\tC*04:01\tC*07:02\t5000\t4900\n",
                encoding="utf-8",
            )
            # 2) arcasHLA JSON for PtA (confirming Class I & adding Class II)
            arc_path = tdir / "PtA.genotype.json"
            arc_path.write_text(
                json.dumps(
                    {
                        "A": ["A*01:01:01", "A*02:01:01"],
                        "B": ["B*27:05:02", "B*44:02:01"],
                        "C": ["C*04:01:01", "C*07:02:01"],
                        "DRB1": ["DRB1*03:01:01", "DRB1*15:01:01"],
                        "DQB1": ["DQB1*02:01:01", "DQB1*06:02:01"],
                    }
                ),
                encoding="utf-8",
            )
            # 3) seq2HLA Class I & II for PtB (homozygous at HLA-A*02:01!)
            seq1_path = tdir / "PtB-ClassI-class.HLAgenotype4digits"
            seq1_path.write_text(
                "#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence\n"
                "A\tA*02:01'\t0.99\tA*02:01\t0.99\n"
                "B\tB*07:02\t0.95\tB*44:02'\t0.94\n"
                "C\tC*07:02\t0.92\tC*05:01\t0.90\n",
                encoding="utf-8",
            )
            # 4) HLA-LA bestguess for PtC
            hlala_path = tdir / "PtC_bestguess_G.txt"
            pd.DataFrame(
                [
                    {"Locus": "A", "Chromosome": 1, "Allele": "HLA-A*03:01:01G", "Q1": 0.99},
                    {"Locus": "A", "Chromosome": 2, "Allele": "HLA-A*24:02:01G", "Q1": 0.99},
                    {"Locus": "B", "Chromosome": 1, "Allele": "HLA-B*07:02:01G", "Q1": 0.98},
                    {"Locus": "B", "Chromosome": 2, "Allele": "HLA-B*35:01:01G", "Q1": 0.97},
                    {"Locus": "C", "Chromosome": 1, "Allele": "HLA-C*04:01:01G", "Q1": 0.99},
                    {"Locus": "C", "Chromosome": 2, "Allele": "HLA-C*04:01:01G", "Q1": 0.99},
                ]
            ).to_csv(hlala_path, sep="\t", index=False)

            res = _run_json(
                "imm_hla_summary.py",
                [
                    "--hla_inputs",
                    str(opt_path),
                    str(arc_path),
                    str(seq1_path),
                    str(hlala_path),
                    "--output_dir",
                    str(tdir / "hla_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "hla-summary-result.json")
            self.assertEqual(res["n_samples"], 3)
            self.assertEqual(res["n_class1_homozygous_samples"], 2)  # PtB (A) and PtC (C)

            hdf = pd.read_csv(res["outputs"]["hla_consensus_genotypes"], sep="\t").set_index("sample_id")
            self.assertEqual(hdf.loc["PtA", "HLA_A_1"], "HLA-A*01:01")
            self.assertEqual(hdf.loc["PtA", "HLA_A_2"], "HLA-A*02:01")
            self.assertEqual(hdf.loc["PtA", "HLA_DRB1_1"], "HLA-DRB1*03:01")
            self.assertFalse(bool(hdf.loc["PtA", "mhc_i_homozygosity_flag"]))
            self.assertGreater(float(hdf.loc["PtA", "mean_class1_hed"]), 0.0)

            self.assertTrue(bool(hdf.loc["PtB", "mhc_i_homozygosity_flag"]))
            self.assertEqual(float(hdf.loc["PtB", "hed_A"]), 0.0)
            self.assertEqual(float(hdf.loc["PtC", "hed_C"]), 0.0)

    def test_5_neoantigen_prioritize_provided_and_pssm(self) -> None:
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # 5a) Upstream MHC binding table
            up_path = tdir / "mhc_preds.tsv"
            pd.DataFrame(
                [
                    {
                        "gene": "TP53",
                        "hla_allele": "HLA-A*02:01",
                        "mt_epitope": "HMTEVVRHC",
                        "wt_epitope": "HMTEVVRRC",
                        "ic50_mt": 28.0,
                        "ic50_wt": 240.0,
                        "tpm": 45.0,
                        "ccf": 0.95,
                    },
                    {
                        "gene": "KRAS",
                        "hla_allele": "HLA-A*02:01",
                        "mt_epitope": "KLVVVGAVGV",
                        "wt_epitope": "KLVVVGAGGV",
                        "ic50_mt": 42.0,
                        "ic50_wt": 3800.0,
                        "tpm": 32.0,
                        "ccf": 0.90,
                    },
                    {
                        "gene": "MUC16_UNEXPR",
                        "hla_allele": "HLA-A*02:01",
                        "mt_epitope": "YLQPRTFLL",
                        "wt_epitope": "YSQPRTFLL",
                        "ic50_mt": 15.0,
                        "ic50_wt": 2500.0,
                        "tpm": 0.1,
                        "ccf": 0.85,
                    },
                    {
                        "gene": "PASSENGER_NONBIND",
                        "hla_allele": "HLA-A*02:01",
                        "mt_epitope": "DRQSTEKPD",
                        "wt_epitope": "DRQSTEKGD",
                        "ic50_mt": 4200.0,
                        "ic50_wt": 4500.0,
                        "tpm": 20.0,
                        "ccf": 0.80,
                    },
                ]
            ).to_csv(up_path, sep="\t", index=False)

            res_up = _run_json(
                "imm_neoantigen_prioritize.py",
                [
                    "--candidates",
                    str(up_path),
                    "--ic50_cutoff_nm",
                    "500",
                    "--min_tpm",
                    "1.0",
                    "--output_dir",
                    str(tdir / "neo_up"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_up, "neoantigen-prioritize-result.json")
            self.assertEqual(res_up["ic50_source"], "provided")
            self.assertEqual(res_up["n_passed_all_gates"], 2)
            self.assertGreaterEqual(res_up["n_strong_binders"], 2)

            ndf = pd.read_csv(res_up["outputs"]["prioritized_neoantigens"], sep="\t").set_index("gene")
            self.assertTrue(bool(ndf.loc["KRAS", "passes_all_gates"]))
            self.assertFalse(bool(ndf.loc["MUC16_UNEXPR", "passes_expression_gate"]))
            self.assertFalse(bool(ndf.loc["PASSENGER_NONBIND", "passes_binding_gate"]))
            self.assertAlmostEqual(float(ndf.loc["KRAS", "agretopicity_dai"]), 3800.0 / 42.0, delta=0.5)
            self.assertEqual(ndf.loc["TP53", "mutation_role"], "tcr_contact")

            # 5b) Somatic flanking peptide table without IC50 -> built-in PSSM estimation
            pep_path = tdir / "somatic_peptides.tsv"
            pd.DataFrame(
                [
                    {
                        "gene": "BRAF_A02",
                        "wt_peptide": "MTEYKLVVVGAGGVGKSALTIQLIQ",
                        "mt_peptide": "MTEYKLVVVGALGVGKSALTIQLIQ",
                        "tpm": 25.0,
                        "vaf": 0.42,
                    },
                    {
                        "gene": "PIK3CA_B07",
                        "wt_peptide": "AISTRDPLSEITEQEKDFLWSHR",
                        "mt_peptide": "AISTRPPLSEITEQEKDFLWSHR",
                        "tpm": 18.0,
                        "vaf": 0.35,
                    },
                ]
            ).to_csv(pep_path, sep="\t", index=False)

            res_pssm = _run_json(
                "imm_neoantigen_prioritize.py",
                [
                    "--candidates",
                    str(pep_path),
                    "--hla_alleles",
                    "HLA-A*02:01",
                    "HLA-B*07:02",
                    "--output_dir",
                    str(tdir / "neo_pssm"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_pssm, "neoantigen-prioritize-result.json")
            self.assertEqual(res_pssm["ic50_source"], "pssm_estimate")
            self.assertEqual(res_pssm["n_candidates"], 2)


if __name__ == "__main__":
    unittest.main()
