"""Contract, offline compilation, and end-to-end runtime tests for the crispr plugin (0.2.0)."""
from __future__ import annotations

import ast
import json
import py_compile
import random
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CRISPR_DIR = REPO_ROOT / "resources" / "plugins" / "crispr"
GD_SKILL = CRISPR_DIR / "skills" / "guide-design-workflow"
SE_SKILL = CRISPR_DIR / "skills" / "screen-editing-workflow"

try:
    import matplotlib  # noqa: F401
    import numpy  # noqa: F401
    import pandas  # noqa: F401
    import scipy  # noqa: F401
    import sklearn  # noqa: F401

    HAS_SCIPY_STACK = True
except ImportError:
    HAS_SCIPY_STACK = False

try:
    import jsonschema  # noqa: F401

    HAS_JSONSCHEMA = True
except ImportError:
    HAS_JSONSCHEMA = False


def _run_json(cmd: list[str]) -> dict:
    proc = subprocess.run([sys.executable, "-B", *cmd], capture_output=True, text=True)
    if proc.returncode != 0:
        raise AssertionError(
            f"Command failed ({proc.returncode}): {' '.join(cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


def _validate_schema_and_sidecars(payload: dict, schema_path: Path) -> None:
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    if HAS_JSONSCHEMA:
        jsonschema.validate(instance=payload, schema=schema)
    else:
        for req_key in schema.get("required", []):
            assert req_key in payload, f"Missing required key '{req_key}' in stdout JSON"
    for art in payload.get("artifacts", []):
        p = Path(art)
        assert p.is_file(), f"Missing artifact file: {art}"
        sidecar = p.with_name(p.name + ".phi-artifact.json")
        assert sidecar.is_file(), f"Missing .phi-artifact.json sidecar for {art}"
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
        assert meta["file"] == p.name
        assert meta["provenance"]["tool"] == payload["tool"]


class TestCrisprPluginContractsAndCompilation(unittest.TestCase):
    def test_layout_manifests_locks_and_no_pycache(self) -> None:
        for reserved in ("mcp", "wrappers", "orchestrator"):
            self.assertFalse((CRISPR_DIR / reserved).exists(), f"Reserved dir {reserved}/ must not exist")

        self.assertEqual(list(CRISPR_DIR.rglob("__pycache__")), [], "No __pycache__ allowed inside plugin")

        manifest = (CRISPR_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("schemaVersion: 1", manifest)
        self.assertIn("id: crispr", manifest)
        self.assertIn("type: plugin", manifest)
        self.assertIn("version: 0.2.0", manifest)
        self.assertIn("toolPrefix: crispr", manifest)
        self.assertIn("- agents/Crispr.md", manifest)
        self.assertIn("- skills/guide-design-workflow", manifest)
        self.assertIn("- skills/screen-editing-workflow", manifest)

        for req_file in ("icon.svg", "README.md", "CHANGELOG.md"):
            self.assertTrue((CRISPR_DIR / req_file).is_file(), f"Missing {req_file}")

        env_yml = (CRISPR_DIR / "environments" / "crispr" / "environment.yml").read_text(encoding="utf-8")
        self.assertIn("name: crispr", env_yml)
        for dep in ("python=3.13", "pandas", "numpy", "scipy", "scikit-learn", "matplotlib-base"):
            self.assertIn(dep, env_yml)

        locks_dir = CRISPR_DIR / "environments" / "crispr" / "locks"
        for plat in ("darwin-arm64.txt", "darwin-x64.txt", "linux-x64.txt"):
            lock_txt = (locks_dir / plat).read_text(encoding="utf-8")
            self.assertIn("# spec: resources/plugins/crispr/environments/crispr/environment.yml", lock_txt)
            non_comments = [l.strip() for l in lock_txt.splitlines() if l.strip() and not l.strip().startswith("#")]
            self.assertEqual(non_comments[0], "@EXPLICIT")
            self.assertGreater(len(non_comments), 10)

    def test_agent_and_skill_contracts(self) -> None:
        agent_txt = (CRISPR_DIR / "agents" / "Crispr.md").read_text(encoding="utf-8")
        self.assertIn("name: Crispr", agent_txt)
        for wt in (
            "wrapper_search",
            "wrapper_inspect",
            "wrapper_plan_run",
            "wrapper_run",
            "wrapper_status",
            "wrapper_wait",
            "wrapper_cancel",
        ):
            self.assertIn(wt, agent_txt)
        self.assertIn("guide-design-workflow", agent_txt)
        self.assertIn("screen-editing-workflow", agent_txt)

        # Check frontmatter even-space sequence indentation and exact script contracts
        skill_specs = [
            (
                GD_SKILL,
                "guide-design-workflow",
                ["design_guides", "offtarget_scan", "library_check", "design_pegrna"],
                [
                    "design-guides-result.json",
                    "offtarget-scan-result.json",
                    "library-check-result.json",
                    "design-pegrna-result.json",
                ],
            ),
            (
                SE_SKILL,
                "screen-editing-workflow",
                ["screen_qc", "screen_hits", "drugz", "editing_quant"],
                [
                    "screen-qc-result.json",
                    "screen-hits-result.json",
                    "drugz-result.json",
                    "editing-quant-result.json",
                ],
            ),
        ]
        for sdir, sname, tool_names, schema_names in skill_specs:
            txt = (sdir / "SKILL.md").read_text(encoding="utf-8")
            end = txt.index("\n---\n", 4)
            fm = txt[4:end]
            self.assertIn(f"name: {sname}", fm)
            self.assertNotRegex(fm, r"^[ ]+toolPrefix:\s*\S+")
            self.assertIn("  environment: plugin:crispr", fm)
            self.assertIn("    - Crispr", fm)
            self.assertTrue((sdir / "references" / "analysis-playbook.md").is_file())
            self.assertIn("references/analysis-playbook.md", txt)

            for line in fm.splitlines():
                m = re.match(r"^( +)-\s", line)
                if m:
                    self.assertEqual(len(m.group(1)) % 2, 0, f"Odd indentation in {sname}: {line}")

            names_found = re.findall(r"^    - name:\s+(\S+)", fm, flags=re.MULTILINE)
            outputs_found = re.findall(r"^      output:\s+(\S+)", fm, flags=re.MULTILINE)
            self.assertEqual(names_found, tool_names)
            self.assertEqual(len(outputs_found), len(tool_names))
            for rel_out in outputs_found:
                self.assertTrue((sdir / rel_out).is_file(), f"Missing output schema {rel_out}")

            for sc in schema_names:
                obj = json.loads((sdir / "schemas" / sc).read_text(encoding="utf-8"))
                self.assertEqual(obj["type"], "object")

    def test_scripts_ast_and_compile(self) -> None:
        scripts = [
            GD_SKILL / "scripts" / "_common.py",
            GD_SKILL / "scripts" / "crispr_design_guides.py",
            GD_SKILL / "scripts" / "crispr_design_pegrna.py",
            GD_SKILL / "scripts" / "crispr_offtarget_scan.py",
            GD_SKILL / "scripts" / "crispr_library_check.py",
            SE_SKILL / "scripts" / "_common.py",
            SE_SKILL / "scripts" / "crispr_screen_qc.py",
            SE_SKILL / "scripts" / "crispr_screen_hits.py",
            SE_SKILL / "scripts" / "crispr_drugz.py",
            SE_SKILL / "scripts" / "crispr_editing_quant.py",
        ]
        with tempfile.TemporaryDirectory() as tmp:
            for p in scripts:
                self.assertTrue(p.is_file(), f"Missing script {p}")
                ast.parse(p.read_text(encoding="utf-8"), filename=str(p))
                py_compile.compile(str(p), cfile=str(Path(tmp) / f"{p.parent.parent.name}_{p.name}.pyc"), doraise=True)


@unittest.skipUnless(HAS_SCIPY_STACK, "Requires numpy/pandas/scipy/scikit-learn/matplotlib")
class TestCrisprPluginRuntime(unittest.TestCase):
    def test_1_design_guides_spcas9_cas12a_and_cbe(self) -> None:
        script = str(GD_SKILL / "scripts" / "crispr_design_guides.py")
        schema = GD_SKILL / "schemas" / "design-guides-result.json"

        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            # Construct a 360 bp synthetic CDS containing:
            # - Early 5' site, middle clean sites, TTTT motif site, BsmBI (CGTCTC) site, Cas12a TTTC site, and CBE CAG codon
            cds = (
                "ATGGCCCTGCCTCTGGAGAGG"  # 21 bp (N-term <5%)
                + "ACTGACCTGGAGATGCTGGTGAGG"  # 24 bp
                + "GCTGATCAGCTGCAGCTGCTGGGG"  # 24 bp (has CAG codon + NGG)
                + "GACGTTTTGCTGACCTGCTGAAGG"  # 24 bp (has TTTT before AGG)
                + "ATCCGTCTCGAGCTGACCTGCTGG"  # 24 bp (has BsmBI CGTCTC before TGG)
                + "TTTCAGCTGATGACCTGCTGAAGTACC"  # 27 bp (has Cas12a TTTC + 23bp)
                + "GCATGCTAGCTGATCGATGCTAGCATCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCTGACTGACTAGCTAG"  # 75 bp
                + "GCTAGCTGACTGATCGTAGCTGACTGATGCTAGCTGACTGATGCTAGCTGACTGATGCTAGCTGACTGATCGTAG"  # 75 bp
                + "CTGATGCTAGCTGACTGATGCTAGCTGACTGATGCTAGCTGACTGATCGATGCTAGCTGACTGA"  # 64 bp
            )
            fa_path = tmp_dir / "target_tp53.fa"
            fa_path.write_text(f">TP53_exon\n{cds}\n", encoding="utf-8")

            # 1a. SpCas9 knockout
            res_ko = _run_json(
                [
                    script,
                    "--target_fasta",
                    str(fa_path),
                    "--nuclease",
                    "SpCas9",
                    "--mode",
                    "knockout",
                    "--cloning_backbone",
                    "lentiCRISPRv2_BsmBI",
                    "--output_dir",
                    str(tmp_dir / "ko_out"),
                ]
            )
            _validate_schema_and_sidecars(res_ko, schema)
            self.assertGreater(res_ko["n_total_candidates"], 10)
            self.assertGreater(res_ko["n_passing_candidates"], 5)

            import pandas as pd

            df_ko = pd.read_csv(res_ko["outputs"]["designed_guides"], sep="\t")
            self.assertTrue((df_ko["has_pol3_terminator"]).any())
            self.assertTrue((df_ko["enzyme_collisions"].fillna("").str.contains("BsmBI")).any())
            self.assertTrue((df_ko["positional_flag"] == "early_nterm_risk").any())
            passing = df_ko[df_ko["passes_filters"]]
            first_row = passing.iloc[0]
            self.assertTrue(first_row["oligo_fwd"].startswith("CACC"))
            self.assertTrue(first_row["oligo_rev"].startswith("AAAC"))
            self.assertTrue(0.0 <= first_row["out_of_frame_score"] <= 100.0)

            # 1b. Cas12a (TTTV 5' PAM, 23 nt spacer)
            res_c12 = _run_json(
                [
                    script,
                    "--target_fasta",
                    str(fa_path),
                    "--nuclease",
                    "Cas12a",
                    "--mode",
                    "knockout",
                    "--output_dir",
                    str(tmp_dir / "c12_out"),
                ]
            )
            _validate_schema_and_sidecars(res_c12, schema)
            df_c12 = pd.read_csv(res_c12["outputs"]["designed_guides"], sep="\t")
            self.assertTrue((df_c12["spacer"].str.len() == 23).all())
            self.assertTrue(df_c12["pam"].str.match(r"^TTT[ACG]$").all())

            # 1c. CBE mode with premature stop codon check (CAG -> TAG)
            # Build frame-0 sequence where CAG sits at pos 4..6 and bystander C sits at pos 8
            # Pos 1..3: ATG, Pos 4..6: CAG (C->T gives TAG stop), Pos 7..9: ACC (bystander C at pos 8), Pos 10..20: TGATGGTGCTG + PAM TGG
            cbe_seq = "ATG" + "CAG" + "ACC" + "TGATGGTGCTG" + "TGG" + "GCTGATGCTGATGCTGATGCTGAT"
            cbe_fa = tmp_dir / "cbe_target.fa"
            cbe_fa.write_text(f">CBE_STOP\n{cbe_seq}\n", encoding="utf-8")
            res_cbe = _run_json(
                [
                    script,
                    "--target_fasta",
                    str(cbe_fa),
                    "--nuclease",
                    "SpCas9",
                    "--mode",
                    "CBE",
                    "--output_dir",
                    str(tmp_dir / "cbe_out"),
                ]
            )
            _validate_schema_and_sidecars(res_cbe, schema)
            df_cbe = pd.read_csv(res_cbe["outputs"]["designed_guides"], sep="\t")
            target_sp = "ATGCAGACCTGATGGTGCTG"
            matched = df_cbe[df_cbe["spacer"] == target_sp]
            self.assertEqual(len(matched), 1)
            self.assertGreaterEqual(int(matched.iloc[0]["target_bases_in_window"]), 2)
            self.assertGreaterEqual(int(matched.iloc[0]["bystander_count"]), 1)
            self.assertTrue(bool(matched.iloc[0]["creates_premature_stop"]))
            self.assertIn("*", str(matched.iloc[0]["predicted_aa_edits"]))

    def test_2_offtarget_scan_mit_cfd_and_nag(self) -> None:
        script = str(GD_SKILL / "scripts" / "crispr_offtarget_scan.py")
        schema = GD_SKILL / "schemas" / "offtarget-scan-result.json"

        g1 = "GACCTGCTGGAGATGCTGGT"  # 20 nt
        g2 = "ATCGTAGCTGATCGATGCTA"  # 20 nt

        # Construct off-targets for g1:
        # 1) Exact match (on-target) with AGG
        # 2) 1mm at pos 1 (5' distal, A vs G) with TGG -> 0 seed mismatches!
        g1_1mm_distal = "AACCTGCTGGAGATGCTGGT"
        # 3) 2mm at pos 2,4 with CGG
        g1_2mm = "GTCTTGCTGGAGATGCTGGT"
        # 4) 1mm at pos 3 with NAG (AAG)
        g1_1mm_nag = "GAGCTGCTGGAGATGCTGGT"

        # Construct sites for g2:
        # 1) Exact on-target + 1 duplicate exact locus -> MIT/CFD <= 50.0
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            guides_tsv = tmp_dir / "guides.tsv"
            guides_tsv.write_text(
                f"guide_id\tspacer\nG1_clean\t{g1}\nG2_dup\t{g2}\n",
                encoding="utf-8",
            )
            pad = "ATATATATATATATATATAT"
            ref_fa = tmp_dir / "ref.fa"
            ref_fa.write_text(
                f">chr1\n{pad}{g1}AGG{pad}{g1_1mm_distal}TGG{pad}{g1_2mm}CGG{pad}{g1_1mm_nag}AAG{pad}\n"
                f">chr2\n{pad}{g2}TGG{pad}{g2}AGG{pad}\n",
                encoding="utf-8",
            )

            res_nag = _run_json(
                [
                    script,
                    "--guides",
                    str(guides_tsv),
                    "--reference_fasta",
                    str(ref_fa),
                    "--max_mismatches",
                    "3",
                    "--include_nag_pam",
                    "--output_dir",
                    str(tmp_dir / "ot_nag"),
                ]
            )
            _validate_schema_and_sidecars(res_nag, schema)
            self.assertTrue(res_nag["include_nag_pam"])
            by_id = {r["guide_id"]: r for r in res_nag["guide_summaries"]}
            self.assertEqual(by_id["G1_clean"]["on_target_hits"], 1)
            self.assertGreaterEqual(by_id["G1_clean"]["n_offtargets"], 3)
            self.assertGreaterEqual(by_id["G1_clean"]["n_zero_seed_mismatch_offtargets"], 1)
            self.assertEqual(by_id["G2_dup"]["on_target_hits"], 2)
            self.assertEqual(by_id["G2_dup"]["n_0mm_dup_offtargets"], 1)
            self.assertLessEqual(by_id["G2_dup"]["mit_specificity_score"], 50.0)
            self.assertLessEqual(by_id["G2_dup"]["cfd_specificity_score"], 50.0)
            self.assertTrue(by_id["G2_dup"]["high_offtarget_risk"])

            # Without NAG PAM (via contract-compliant --canonical_pam_only flag) -> fewer off-targets for G1_clean
            res_no_nag = _run_json(
                [
                    script,
                    "--guides",
                    str(guides_tsv),
                    "--reference_fasta",
                    str(ref_fa),
                    "--max_mismatches",
                    "3",
                    "--canonical_pam_only",
                    "--output_dir",
                    str(tmp_dir / "ot_nonag"),
                ]
            )
            _validate_schema_and_sidecars(res_no_nag, schema)
            self.assertFalse(res_no_nag["include_nag_pam"])
            by_id2 = {r["guide_id"]: r for r in res_no_nag["guide_summaries"]}
            self.assertLess(by_id2["G1_clean"]["n_offtargets"], by_id["G1_clean"]["n_offtargets"])

            # Cas12a 5' PAM (TTTV) off-target scan
            c12_sp = "AGCTGATGACCTGCTGAAGTACC"  # 23 nt
            c12_ot = "AGCTGATGACCTGCTGAAGTAGC"  # 1 mismatch at 3' PAM-distal pos 22 -> 0 5'-seed mismatches
            c12_guides = tmp_dir / "c12_guides.tsv"
            c12_guides.write_text(f"guide_id\tspacer\tnuclease\nC12_g1\t{c12_sp}\tCas12a\n", encoding="utf-8")
            c12_ref = tmp_dir / "c12_ref.fa"
            c12_ref.write_text(f">chrC12\n{pad}TTTC{c12_sp}{pad}TTTA{c12_ot}{pad}\n", encoding="utf-8")
            res_c12_ot = _run_json(
                [
                    script,
                    "--guides",
                    str(c12_guides),
                    "--reference_fasta",
                    str(c12_ref),
                    "--max_mismatches",
                    "2",
                    "--output_dir",
                    str(tmp_dir / "ot_c12"),
                ]
            )
            _validate_schema_and_sidecars(res_c12_ot, schema)
            c12_sum = res_c12_ot["guide_summaries"][0]
            self.assertEqual(c12_sum["on_target_hits"], 1)
            self.assertEqual(c12_sum["n_1mm_offtargets"], 1)
            self.assertEqual(c12_sum["n_zero_seed_mismatch_offtargets"], 1)

    def test_3_library_check_catches_defects(self) -> None:
        script = str(GD_SKILL / "scripts" / "crispr_library_check.py")
        schema = GD_SKILL / "schemas" / "library-check-result.json"

        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            lib_path = tmp_dir / "lib.tsv"
            rows = [
                "sgRNA\tseq\tgene",
                "sg1\tGACCTGCTGGAGATGCTGGT\tTP53",
                "sg1\tATGACCAGCTGATGGTGCTG\tTP53",  # duplicate ID
                "sg3\tGACCTGCTGGAGATGCTGGT\tTP53",  # duplicate sequence
                "sg4\tGACGTTTTGCTGACCTGCTG\tKRAS",  # TTTT Pol III terminator + <3 guides for KRAS
                "sg5\tATCCGTCTCGAGCTGACCTG\tKRAS",  # BsmBI CGTCTC collision
                "sg6\tATCCGTCTCGAGCTGACNTG\tBRAF",  # invalid base N + <3 guides for BRAF
                "ntc1\tGCTAGCTGATCGATGCTAGC\tNonTargetingControl",
            ]
            lib_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

            res = _run_json(
                [
                    script,
                    "--library",
                    str(lib_path),
                    "--output_dir",
                    str(tmp_dir / "lib_qc"),
                ]
            )
            _validate_schema_and_sidecars(res, schema)
            self.assertEqual(res["qc_status"], "fail")
            self.assertEqual(res["n_duplicate_ids"], 2)
            self.assertEqual(res["n_duplicate_sequences"], 2)
            self.assertEqual(res["n_invalid_sequences"], 1)
            self.assertGreaterEqual(res["n_pol3_tttt"], 1)
            self.assertGreaterEqual(res["n_enzyme_collisions"], 2)
            self.assertEqual(res["n_control_guides"], 1)
            self.assertGreaterEqual(res["n_genes_low_coverage"], 2)

    def test_4_screen_qc_gini_and_nnmd(self) -> None:
        script = str(SE_SKILL / "scripts" / "crispr_screen_qc.py")
        schema = SE_SKILL / "schemas" / "screen-qc-result.json"

        rng = random.Random(101)
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            ess_genes = [f"ESS_{i}" for i in range(12)]
            non_genes = [f"NON_{i}" for i in range(12)]
            bg_genes = [f"GENE_{i}" for i in range(50)]

            lines = ["sgRNA\tGene\tT0_Plasmid\tDay14_Rep1\tDay14_Rep2\tBottlenecked_Sample"]
            for g in ess_genes + non_genes + bg_genes:
                for k in range(4):
                    sg = f"{g}_sg{k+1}"
                    t0 = rng.randint(450, 580)
                    if g.startswith("ESS_"):
                        r1 = rng.randint(35, 75)
                        r2 = rng.randint(40, 80)
                    else:
                        r1 = rng.randint(420, 560)
                        r2 = rng.randint(410, 570)
                    # Extreme skew in Bottlenecked_Sample (90% zeros, huge spikes on few sgRNAs)
                    bot = 8000 if (k == 0 and g.startswith("GENE_0")) else (0 if rng.random() < 0.7 else 15)
                    lines.append(f"{sg}\t{g}\t{t0}\t{r1}\t{r2}\t{bot}")

            counts_path = tmp_dir / "counts.txt"
            counts_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            ess_path = tmp_dir / "essential.txt"
            non_path = tmp_dir / "nonessential.txt"
            ess_path.write_text("\n".join(ess_genes) + "\n", encoding="utf-8")
            non_path.write_text("\n".join(non_genes) + "\n", encoding="utf-8")

            res = _run_json(
                [
                    script,
                    "--counts",
                    str(counts_path),
                    "--essential_genes",
                    str(ess_path),
                    "--nonessential_genes",
                    str(non_path),
                    "--output_dir",
                    str(tmp_dir / "qc_out"),
                ]
            )
            _validate_schema_and_sidecars(res, schema)
            self.assertEqual(res["reference_sample"], "T0_Plasmid")
            by_s = {s["sample"]: s for s in res["samples"]}
            self.assertLess(by_s["T0_Plasmid"]["gini_index"], 0.15)
            self.assertGreater(by_s["Bottlenecked_Sample"]["gini_index"], 0.30)
            self.assertIn("severe_bottlenecking", by_s["Bottlenecked_Sample"]["qc_flags"])

            nnmd_by_s = {n["sample"]: n for n in res["nnmd_summaries"]}
            self.assertLess(nnmd_by_s["Day14_Rep1"]["nnmd"], -1.25)
            self.assertGreater(nnmd_by_s["Day14_Rep1"]["essentiality_roc_auc"], 0.90)
            self.assertTrue(nnmd_by_s["Day14_Rep1"]["strong_separation"])

    def test_5_screen_hits_recovers_depleted_and_enriched_genes(self) -> None:
        script = str(SE_SKILL / "scripts" / "crispr_screen_hits.py")
        schema = SE_SKILL / "schemas" / "screen-hits-result.json"

        rng = random.Random(2026)
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            dep_truth = {"MYC", "RPL5", "PSMA1"}
            enr_truth = {"TP53", "NF2"}
            genes = list(dep_truth) + list(enr_truth) + [f"NEUT_{i:03d}" for i in range(110)] + ["NTC"]

            lines = ["sgRNA\tGene\tCtrl_1\tCtrl_2\tTrt_1\tTrt_2"]
            for g in genes:
                n_sg = 10 if g == "NTC" else 4
                for k in range(n_sg):
                    c1 = rng.randint(420, 620)
                    c2 = rng.randint(420, 620)
                    if g in dep_truth:
                        t1 = rng.randint(25, 65)
                        t2 = rng.randint(25, 70)
                    elif g in enr_truth:
                        t1 = rng.randint(2400, 3200)
                        t2 = rng.randint(2300, 3100)
                    else:
                        t1 = max(20, int(c1 * rng.uniform(0.85, 1.15)))
                        t2 = max(20, int(c2 * rng.uniform(0.85, 1.15)))
                    lines.append(f"{g}_sg{k+1}\t{g}\t{c1}\t{c2}\t{t1}\t{t2}")

            counts_path = tmp_dir / "screen_counts.tsv"
            counts_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

            res = _run_json(
                [
                    script,
                    "--counts",
                    str(counts_path),
                    "--treatment_cols",
                    "Trt_1",
                    "--treatment_cols",
                    "Trt_2",
                    "--control_cols",
                    "Ctrl_1",
                    "--control_cols",
                    "Ctrl_2",
                    "--norm_method",
                    "control_sgrna",
                    "--n_permutations",
                    "1000",
                    "--output_dir",
                    str(tmp_dir / "hits_out"),
                ]
            )
            _validate_schema_and_sidecars(res, schema)
            self.assertEqual(res["norm_method"], "control_sgrna")
            dep_hits = {r["gene"] for r in res["top_depleted_genes"] if r["neg_fdr"] < 0.05}
            enr_hits = {r["gene"] for r in res["top_enriched_genes"] if r["pos_fdr"] < 0.05}
            self.assertTrue(dep_truth.issubset(dep_hits), f"Expected {dep_truth} in {dep_hits}")
            self.assertTrue(enr_truth.issubset(enr_hits), f"Expected {enr_truth} in {enr_hits}")

    def test_6_editing_quant_alleles_cbe_and_crispresso2(self) -> None:
        script = str(SE_SKILL / "scripts" / "crispr_editing_quant.py")
        schema = SE_SKILL / "schemas" / "editing-quant-result.json"

        # 60 bp WT amplicon with SpCas9 guide at indices 18..38 (20 nt) + PAM TGG at 38..41, cut at 35
        # Guide has C at pos 6 (idx 23) and C at pos 8 (idx 25) for CBE testing too!
        prefix = "GCTAGCTGATCGATGCTA"  # 18 bp
        guide = "ATGACCACCTGATGGTGCTG"  # 20 nt: pos 6=C, pos 7=C, pos 8=A; let's make pos 6=C, pos 8=C:
        guide = "ATGACCTCCTGATGGTGCTG"  # pos1=A,2=T,3=G,4=A,5=C,6=C,7=T,8=C...
        pam_tail = "TGGGCTAGCTGATCGATGCTAG"  # 22 bp
        wt_amp = prefix + guide + pam_tail  # 60 bp
        # HDR template: introduces 2 silent SNVs near cut site (indices 34, 37)
        hdr_amp = wt_amp[:34] + "AT" + wt_amp[36:]

        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)

            # 6a. Mode A: Nuclease KO + HDR
            fs_del2 = wt_amp[:33] + wt_amp[35:]  # -2 bp frameshift deletion around cut 35
            inf_del3 = wt_amp[:32] + wt_amp[35:]  # -3 bp in-frame deletion
            imp_hdr = hdr_amp[:32] + hdr_amp[33:]  # HDR SNVs + extra -1 bp deletion

            alleles_a = tmp_dir / "alleles_nuclease.tsv"
            alleles_a.write_text(
                "\n".join(
                    [
                        "sample\tsequence\tcount",
                        f"EditSample\t{wt_amp}\t400",
                        f"EditSample\t{fs_del2}\t300",
                        f"EditSample\t{inf_del3}\t100",
                        f"EditSample\t{hdr_amp}\t150",
                        f"EditSample\t{imp_hdr}\t50",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            res_a = _run_json(
                [
                    script,
                    "--alleles",
                    str(alleles_a),
                    "--amplicon_seq",
                    wt_amp,
                    "--guide_seq",
                    guide,
                    "--hdr_seq",
                    hdr_amp,
                    "--editor_mode",
                    "nuclease",
                    "--output_dir",
                    str(tmp_dir / "eq_nuclease"),
                ]
            )
            _validate_schema_and_sidecars(res_a, schema)
            s_a = res_a["samples"][0]
            self.assertEqual(s_a["total_reads"], 1000)
            self.assertAlmostEqual(s_a["unmodified_pct"], 40.0, places=2)
            self.assertAlmostEqual(s_a["frameshift_indel_pct"], 30.0, places=2)
            self.assertAlmostEqual(s_a["in_frame_indel_pct"], 10.0, places=2)
            self.assertAlmostEqual(s_a["hdr_pct"], 15.0, places=2)
            self.assertAlmostEqual(s_a["imperfect_hdr_pct"], 5.0, places=2)

            # 6b. Mode A: CBE base editing (C6->T target vs C8->T bystander + indel byproduct)
            # Guide indices inside wt_amp: g_start=18; pos6 = 18+5=23 ('C'), pos8 = 18+7=25 ('C')
            cbe_pure = wt_amp[:23] + "T" + wt_amp[24:]
            cbe_byst = wt_amp[:23] + "T" + wt_amp[24:25] + "T" + wt_amp[26:]
            alleles_cbe = tmp_dir / "alleles_cbe.tsv"
            alleles_cbe.write_text(
                "\n".join(
                    [
                        "sample\tsequence\tcount",
                        f"CBE_Run\t{wt_amp}\t200",
                        f"CBE_Run\t{cbe_pure}\t600",
                        f"CBE_Run\t{cbe_byst}\t150",
                        f"CBE_Run\t{wt_amp[:34] + wt_amp[35:]}\t50",  # 1 bp deletion byproduct
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            res_cbe = _run_json(
                [
                    script,
                    "--alleles",
                    str(alleles_cbe),
                    "--amplicon_seq",
                    wt_amp,
                    "--guide_seq",
                    guide,
                    "--editor_mode",
                    "CBE",
                    "--window_bp",
                    "15",
                    "--output_dir",
                    str(tmp_dir / "eq_cbe"),
                ]
            )
            _validate_schema_and_sidecars(res_cbe, schema)
            s_cbe = res_cbe["samples"][0]
            self.assertAlmostEqual(s_cbe["max_be_target_conversion_pct"], 75.0, places=2)
            self.assertAlmostEqual(s_cbe["be_bystander_edit_pct"], 15.0, places=2)
            self.assertAlmostEqual(s_cbe["indel_byproduct_pct"], 5.0, places=2)

            # 6c. Mode B: CRISPResso2 Alleles_frequency_table_around_sgRNA_*.txt
            sub_ref = wt_amp[15:48]  # 33 bp sub-window covering guide & cut site
            sub_fs4 = sub_ref[:16] + "----" + sub_ref[20:]
            sub_inf3 = sub_ref[:19] + "AAA" + sub_ref[19:]
            sub_ref_ins = sub_ref[:19] + "---" + sub_ref[19:]
            crispresso_tsv = tmp_dir / "Alleles_frequency_table_around_sgRNA_demo.txt"
            crispresso_tsv.write_text(
                "\n".join(
                    [
                        "Aligned_Sequence\tReference_Sequence\tRead_Status\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads",
                        f"{sub_ref}\t{sub_ref}\tUNMODIFIED\t0\t0\t0\t500\t50.0",
                        f"{sub_fs4}\t{sub_ref}\tMODIFIED\t4\t0\t0\t350\t35.0",
                        f"{sub_inf3}\t{sub_ref_ins}\tMODIFIED\t0\t3\t0\t150\t15.0",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            res_b = _run_json(
                [
                    script,
                    "--crispresso_table",
                    str(crispresso_tsv),
                    "--amplicon_seq",
                    wt_amp,
                    "--guide_seq",
                    guide,
                    "--output_dir",
                    str(tmp_dir / "eq_crispresso"),
                ]
            )
            _validate_schema_and_sidecars(res_b, schema)
            self.assertEqual(res_b["input_mode"], "crispresso2")
            s_b = res_b["samples"][0]
            self.assertAlmostEqual(s_b["unmodified_pct"], 50.0, places=2)
            self.assertAlmostEqual(s_b["frameshift_indel_pct"], 35.0, places=2)
            self.assertAlmostEqual(s_b["in_frame_indel_pct"], 15.0, places=2)

    def test_7_design_pegrna_pe3b_and_epgrna(self) -> None:
        script = str(GD_SKILL / "scripts" / "crispr_design_pegrna.py")
        schema = GD_SKILL / "schemas" / "design-pegrna-result.json"

        # 120 bp WT CDS with SpCas9 PAM (AGG) right upstream of edit + '-' strand CCN enabling PE3b
        left_flank = "ATGGCTGACCTGCTGGAGATGCTGGTGAAGCTGATCGATGCTAGCGACCTGCTGGAG"  # 57 bp
        # Spacer region: 20 nt + PAM AGG; edit changes G->T 2 bp downstream of nick (disrupting seed!)
        # Also creates a '-' strand PAM / seed match for PE3b!
        core_wt = "ATGCAGACCTGATGGTGCTGAGG" + "CCTAGCTGACCTGCTGGAGATGCTGGTGAAGCTGAT"
        core_ed = "ATGCAGACCTGATGGTGCATAGG" + "CCTAGCTGACCTGCTGGAGATGCTGGTGAAGCTGAT"
        wt_seq = left_flank + core_wt
        ed_seq = left_flank + core_ed

        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            res = _run_json(
                [
                    script,
                    "--wt_seq",
                    wt_seq,
                    "--edited_seq",
                    ed_seq,
                    "--add_epgrna_motif",
                    "--cds_frame_offset",
                    "0",
                    "--output_dir",
                    str(tmp_dir / "pegrna_out"),
                ]
            )
            _validate_schema_and_sidecars(res, schema)
            self.assertEqual(res["edit_type"], "substitution")
            self.assertTrue(res["epeg_motif_added"])
            self.assertGreater(res["n_pegrna_candidates"], 5)
            self.assertGreater(res["n_ngrna_candidates"], 0)
            top_peg = res["top_pegrnas"][0]
            self.assertIn("CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA", top_peg["extension_3prime"])

    def test_8_drugz_chemogenetic_screen(self) -> None:
        script = str(SE_SKILL / "scripts" / "crispr_drugz.py")
        schema = SE_SKILL / "schemas" / "drugz-result.json"

        rng = random.Random(2027)
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            sl_truth = {"BRCA1", "BRCA2", "PALB2"}
            sup_truth = {"TP53BP1", "SHLD2"}
            neut_genes = [f"GENE_{i:03d}" for i in range(120)]
            lines = ["sgRNA\tGene\tDMSO_1\tDMSO_2\tPARPi_1\tPARPi_2"]

            for g in list(sl_truth) + list(sup_truth) + neut_genes + ["PAN_ESS"]:
                for k in range(4):
                    base = rng.randint(250, 900)
                    d1 = max(20, int(base * rng.uniform(0.92, 1.08)))
                    d2 = max(20, int(base * rng.uniform(0.92, 1.08)))
                    if g in sl_truth:
                        p1 = max(5, int(d1 * rng.uniform(0.12, 0.22)))
                        p2 = max(5, int(d2 * rng.uniform(0.12, 0.22)))
                    elif g in sup_truth:
                        p1 = int(d1 * rng.uniform(4.2, 5.5))
                        p2 = int(d2 * rng.uniform(4.2, 5.5))
                    elif g == "PAN_ESS":
                        d1, d2, p1, p2 = 12, 10, 11, 9
                    else:
                        p1 = max(20, int(d1 * rng.uniform(0.88, 1.12)))
                        p2 = max(20, int(d2 * rng.uniform(0.88, 1.12)))
                    lines.append(f"{g}_sg{k+1}\t{g}\t{d1}\t{d2}\t{p1}\t{p2}")

            counts_path = tmp_dir / "parpi_counts.tsv"
            counts_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

            res = _run_json(
                [
                    script,
                    "--counts",
                    str(counts_path),
                    "--drug_cols",
                    "PARPi_1",
                    "--drug_cols",
                    "PARPi_2",
                    "--control_cols",
                    "DMSO_1",
                    "--control_cols",
                    "DMSO_2",
                    "--remove_genes",
                    "PAN_ESS",
                    "--output_dir",
                    str(tmp_dir / "drugz_out"),
                ]
            )
            _validate_schema_and_sidecars(res, schema)
            self.assertEqual(res["mode"], "paired")
            sl_hits = {r["gene"] for r in res["top_synthetic_lethal"] if r["fdr_synth"] < 0.05}
            sup_hits = {r["gene"] for r in res["top_suppressors"] if r["fdr_supp"] < 0.05}
            self.assertTrue(sl_truth.issubset(sl_hits), f"Expected {sl_truth} in {sl_hits}")
            self.assertTrue(sup_truth.issubset(sup_hits), f"Expected {sup_truth} in {sup_hits}")


if __name__ == "__main__":
    unittest.main()
