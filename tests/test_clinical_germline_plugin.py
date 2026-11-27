#!/usr/bin/env python3
"""Contract, layout, AST, schema, and end-to-end runtime tests for clinical-germline (0.1.0)."""

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
PLUGIN_DIR = ROOT / "resources" / "plugins" / "clinical-germline"
SKILL_DIR = PLUGIN_DIR / "skills" / "germline-workflow"

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
        testcase.assertEqual(meta["provenance"]["skill"], "germline-workflow")


class TestClinicalGermlineContract(unittest.TestCase):
    def test_plugin_layout_locks_schemas_and_ast(self) -> None:
        manifest = (PLUGIN_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("schemaVersion: 1", manifest)
        self.assertIn("id: clinical-germline", manifest)
        self.assertIn("type: plugin", manifest)
        self.assertIn("version: 0.1.0", manifest)
        self.assertIn("toolPrefix: germ", manifest)
        self.assertIn("- agents/ClinicalGermline.md", manifest)
        self.assertIn("- skills/germline-workflow", manifest)
        self.assertIn("spec: environments/germ/environment.yml", manifest)

        m_sum = re.search(r"^summary:\s*>-\n((?:  .*\n)+)", manifest, re.M)
        self.assertIsNotNone(m_sum)
        summary_text = " ".join(line.strip() for line in m_sum.group(1).splitlines())
        self.assertLessEqual(len(summary_text), 260, f"summary too long: {len(summary_text)}")

        for req in ("icon.svg", "README.md", "CHANGELOG.md", "agents/ClinicalGermline.md"):
            self.assertTrue((PLUGIN_DIR / req).is_file(), f"missing {req}")
        for reserved in ("mcp", "wrappers", "orchestrator"):
            self.assertFalse((PLUGIN_DIR / reserved).exists(), f"reserved dir {reserved} exists")

        readme_txt = (PLUGIN_DIR / "README.md").read_text(encoding="utf-8")
        for ref in ("ensemblvep", "snpeff", "clair3", "deepvariant", "manta", "nf-core/sarek", "ols4"):
            self.assertIn(ref, readme_txt)

        agent_txt = (PLUGIN_DIR / "agents/ClinicalGermline.md").read_text(encoding="utf-8")
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
            lock = PLUGIN_DIR / f"environments/germ/locks/{plat}.txt"
            self.assertTrue(lock.is_file(), f"missing lock {plat}")
            ltxt = lock.read_text(encoding="utf-8")
            self.assertIn("@EXPLICIT", ltxt)
            self.assertIn("# spec: resources/plugins/clinical-germline/environments/germ/environment.yml", ltxt)

        stext = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
        sfm = stext[4 : stext.index("\n---\n", 4)]
        sbody = stext[stext.index("\n---\n", 4) + 5 :]
        self.assertIn("name: germline-workflow", sfm)
        self.assertIn("  environment: plugin:germ", sfm)
        self.assertIn("    - ClinicalGermline", sfm)
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
                "vcf_filter",
                "acmg_prioritize",
                "phenotype_rank",
                "pgx_star_allele",
                "cnv_sv_integrate",
            ],
        )
        self.assertEqual(len(outs), 5)
        for rel_out in outs:
            schema_file = SKILL_DIR / rel_out
            self.assertTrue(schema_file.is_file(), f"missing schema {rel_out}")
            json.loads(schema_file.read_text(encoding="utf-8"))

        exp_scripts = {
            "_common.py",
            "germ_vcf_filter.py",
            "germ_acmg_prioritize.py",
            "germ_phenotype_rank.py",
            "germ_pgx_star_allele.py",
            "germ_cnv_sv_integrate.py",
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
class TestClinicalGermlineRuntimeE2E(unittest.TestCase):
    def test_1_vcf_filter_trio_and_singleton(self) -> None:
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            vcf_path = tdir / "trio_annotated.vcf"
            vcf_lines = [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|HGVSc|HGVSp|gnomAD_AF|CLIN_SIG|REVEL|CADD_PHRED|SpliceAI_max|pLI|LOEUF">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND\tMOTHER\tFATHER",
                # 1) SCN1A de novo stop_gained
                "chr2\t166053000\t.\tG\tA\t999\tPASS\tCSQ=A|stop_gained|HIGH|SCN1A|c.1234C>T|p.Arg412*|0.0|pathogenic||38.0|0.02|1.0|0.12\tGT:AD:DP:GQ\t0/1:24,22:46:99\t0/0:40,0:40:99\t0/0:44,0:44:99",
                # 2) ATP7B phased in trans compound het (maternal + paternal)
                "chr13\t51958000\t.\tC\tT\t800\tPASS\tCSQ=T|missense_variant|MODERATE|ATP7B|c.2333G>A|p.Arg778Leu|0.0004|pathogenic|0.91|27.5|0.01|0.02|0.65\tGT:AD:DP:GQ\t0/1:20,21:41:99\t0/1:22,18:40:99\t0/0:38,0:38:99",
                "chr13\t51975000\t.\tG\tC\t820\tPASS\tCSQ=C|missense_variant|MODERATE|ATP7B|c.2975C>G|p.Pro992Arg|0.0002|pathogenic|0.88|26.4|0.01|0.02|0.65\tGT:AD:DP:GQ\t0/1:19,20:39:99\t0/0:42,0:42:99\t0/1:21,19:40:99",
                # 3) TTN in-cis pair (both maternal only -> excluded from compound_heterozygous)
                "chr2\t178600100\t.\tA\tG\t500\tPASS\tCSQ=G|missense_variant|MODERATE|TTN|c.1000T>C|p.Tyr334His|0.002||0.35|18.0|0.01|0.0|0.9\tGT:AD:DP:GQ\t0/1:20,18:38:90\t0/1:18,17:35:90\t0/0:36,0:36:95",
                "chr2\t178610200\t.\tC\tT\t500\tPASS\tCSQ=T|missense_variant|MODERATE|TTN|c.2000G>A|p.Arg667Gln|0.003||0.31|16.5|0.01|0.0|0.9\tGT:AD:DP:GQ\t0/1:22,19:41:92\t0/1:19,20:39:91\t0/0:40,0:40:96",
                # 4) GAA autosomal recessive homozygous
                "chr17\t80105000\t.\tG\tA\t900\tPASS\tCSQ=A|missense_variant|MODERATE|GAA|c.546G>A|p.Thr182Thr|0.001|pathogenic|0.79|26.0|0.65|0.0|0.8\tGT:AD:DP:GQ\t1/1:0,42:42:99\t0/1:20,19:39:99\t0/1:21,22:43:99",
                # 5) BRCA1 rare dominant / ACMG SF v3.2
                "chr17\t43091000\t.\tCT\tC\t950\tPASS\tCSQ=C|frameshift_variant|HIGH|BRCA1|c.5266dup|p.Gln1756Profs|0.00002|pathogenic||35.0|0.0|0.99|0.18\tGT:AD:DP:GQ\t0/1:23,21:44:99\t0/1:25,20:45:99\t0/0:41,0:41:99",
                # 6) Noise: low DP & synonymous
                "chr1\t100000\t.\tA\tG\t50\tPASS\tCSQ=G|synonymous_variant|LOW|NOC2L|c.30A>G|p.Leu10=|0.00001||0.10|5.0|0.0|0.1|0.8\tGT:AD:DP:GQ\t0/1:3,2:5:12\t0/0:20,0:20:80\t0/0:20,0:20:80",
            ]
            vcf_path.write_text("\n".join(vcf_lines) + "\n", encoding="utf-8")

            res_trio = _run_json(
                "germ_vcf_filter.py",
                [
                    "--variants", str(vcf_path),
                    "--proband_id", "PROBAND",
                    "--mother_id", "MOTHER",
                    "--father_id", "FATHER",
                    "--proband_sex", "female",
                    "--output_dir", str(tdir / "trio_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_trio, "vcf-filter-result.json")
            self.assertTrue(res_trio["trio_mode"])
            self.assertEqual(res_trio["total_variants_input"], 8)

            cdf = pd.read_csv(res_trio["outputs"]["filtered_germline_candidates"], sep="\t")
            scn1a = cdf[cdf["gene"] == "SCN1A"].iloc[0]
            self.assertEqual(scn1a["inheritance_mode"], "de_novo")
            self.assertEqual(scn1a["phasing_status"], "trio_de_novo")

            atp7b = cdf[cdf["gene"] == "ATP7B"]
            self.assertEqual(len(atp7b), 2)
            self.assertTrue((atp7b["inheritance_mode"] == "compound_heterozygous").all())
            self.assertTrue((atp7b["phasing_status"] == "phased_in_trans").all())

            # TTN in-cis pair must NOT be called compound_heterozygous
            self.assertFalse((cdf[cdf["gene"] == "TTN"]["inheritance_mode"] == "compound_heterozygous").any())
            self.assertEqual(cdf[cdf["gene"] == "GAA"].iloc[0]["inheritance_mode"], "autosomal_recessive_hom")

            # Singleton mode test: ATP7B should be flagged unphased_candidate_comp_het
            res_single = _run_json(
                "germ_vcf_filter.py",
                [
                    "--variants", str(vcf_path),
                    "--proband_id", "PROBAND",
                    "--output_dir", str(tdir / "single_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_single, "vcf-filter-result.json")
            self.assertFalse(res_single["trio_mode"])
            sdf = pd.read_csv(res_single["outputs"]["filtered_germline_candidates"], sep="\t")
            atp_s = sdf[sdf["gene"] == "ATP7B"]
            self.assertTrue((atp_s["inheritance_mode"] == "unphased_candidate_comp_het").all())
            self.assertTrue((atp_s["phasing_status"] == "unphased_singleton").all())

    def test_2_and_3_acmg_prioritize_and_phenotype_rank(self) -> None:
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            cand_path = tdir / "candidates.tsv"
            pd.DataFrame(
                [
                    {
                        "variant_id": "chr2:166053000:G>A",
                        "gene": "SCN1A",
                        "consequence": "stop_gained",
                        "gnomad_af": 0.0,
                        "clinvar_sig": "",
                        "revel": None,
                        "cadd_phred": 38.0,
                        "spliceai_max": 0.02,
                        "pli": 1.0,
                        "loeuf": 0.12,
                        "inheritance_mode": "de_novo",
                        "phasing_status": "trio_de_novo",
                    },
                    {
                        "variant_id": "chr17:43091000:CT>C",
                        "gene": "BRCA1",
                        "consequence": "frameshift_variant",
                        "gnomad_af": 0.000005,
                        "clinvar_sig": "Pathogenic",
                        "revel": None,
                        "cadd_phred": 35.0,
                        "spliceai_max": 0.0,
                        "pli": 0.99,
                        "loeuf": 0.18,
                        "inheritance_mode": "autosomal_dominant",
                        "phasing_status": "maternal_only",
                    },
                    {
                        "variant_id": "chr13:51958000:C>T",
                        "gene": "ATP7B",
                        "consequence": "missense_variant",
                        "gnomad_af": 0.00005,
                        "clinvar_sig": "",
                        "revel": 0.85,
                        "cadd_phred": 27.0,
                        "spliceai_max": 0.01,
                        "pli": 0.02,
                        "loeuf": 0.65,
                        "inheritance_mode": "compound_heterozygous",
                        "phasing_status": "phased_in_trans",
                    },
                    {
                        "variant_id": "chr11:47340000:G>A",
                        "gene": "MYBPC3",
                        "consequence": "missense_variant",
                        "gnomad_af": 0.0002,
                        "clinvar_sig": "Uncertain_significance",
                        "revel": 0.45,
                        "cadd_phred": 19.0,
                        "spliceai_max": 0.01,
                        "pli": 0.1,
                        "loeuf": 0.6,
                        "inheritance_mode": "autosomal_dominant",
                        "phasing_status": "unphased_singleton",
                    },
                    {
                        "variant_id": "chr1:900000:A>G",
                        "gene": "OR4F5",
                        "consequence": "synonymous_variant",
                        "gnomad_af": 0.12,
                        "clinvar_sig": "Benign",
                        "revel": 0.10,
                        "cadd_phred": 4.0,
                        "spliceai_max": 0.0,
                        "pli": 0.0,
                        "loeuf": 1.2,
                        "inheritance_mode": "autosomal_dominant",
                        "phasing_status": "unphased_singleton",
                    },
                ]
            ).to_csv(cand_path, sep="\t", index=False)

            res_acmg = _run_json(
                "germ_acmg_prioritize.py",
                ["--variants", str(cand_path), "--output_dir", str(tdir / "acmg_out")],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_acmg, "acmg-prioritize-result.json")
            self.assertGreaterEqual(res_acmg["n_pathogenic_or_likely"], 2)
            self.assertEqual(res_acmg["n_acmg_sf_reportable"], 1)

            adf = pd.read_csv(res_acmg["outputs"]["acmg_classified_variants"], sep="\t").set_index("gene")
            # SCN1A: PVS1(+8) + PS2(+4) + PM2_Supporting(+1) = 13 -> Pathogenic
            self.assertEqual(int(adf.loc["SCN1A", "acmg_points"]), 13)
            self.assertEqual(adf.loc["SCN1A", "acmg_class"], "Pathogenic")
            self.assertGreater(float(adf.loc["SCN1A", "tavtigian_posterior"]), 0.99)
            self.assertIn("PVS1", str(adf.loc["SCN1A", "acmg_evidence"]))
            self.assertIn("PS2", str(adf.loc["SCN1A", "acmg_evidence"]))

            self.assertTrue(bool(adf.loc["BRCA1", "acmg_sf_reportable"]))
            self.assertIn("PM3", str(adf.loc["ATP7B", "acmg_evidence"]))
            self.assertEqual(adf.loc["OR4F5", "acmg_class"], "Benign")

            # Step 3: HPO Phenotype Ranking for Epilepsy (SCN1A must rank #1!)
            res_pheno = _run_json(
                "germ_phenotype_rank.py",
                [
                    "--variants", res_acmg["outputs"]["acmg_classified_variants"],
                    "--hpo_terms", "HP:0001250", "HP:0002373", "HP:0001263",
                    "--output_dir", str(tdir / "pheno_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_pheno, "phenotype-rank-result.json")
            gdf = pd.read_csv(res_pheno["outputs"]["diagnostic_gene_ranking"], sep="\t")
            self.assertEqual(gdf.iloc[0]["gene"], "SCN1A")
            self.assertGreater(float(gdf.iloc[0]["phenotype_sim"]), 0.70)
            self.assertGreater(float(gdf.iloc[0]["combined_diagnostic_score"]), 85.0)

    def test_4_pgx_star_allele_cpic_alerts(self) -> None:
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            pgx_path = tdir / "patient_pgx_rsids.tsv"
            pd.DataFrame(
                [
                    {"sample_id": "PT01", "rsid": "rs4244285", "genotype": "A/A"},   # CYP2C19 *2/*2 -> PM -> Clopidogrel alternative_therapy
                    {"sample_id": "PT01", "rsid": "rs3918290", "genotype": "G/A"},   # DPYD *1/*2A -> IM (AS=1.0) -> 5-FU dose_reduction_required
                    {"sample_id": "PT01", "rsid": "rs3892097", "genotype": "T/T"},   # CYP2D6 *4/*4 -> PM -> Codeine contraindicated_or_avoid
                    {"sample_id": "PT01", "rsid": "rs9923231", "genotype": "G/A"},   # VKORC1 Decreased_Function -> Warfarin dose_reduction_required
                    {"sample_id": "PT01", "rsid": "rs4149056", "genotype": "T/C"},   # SLCO1B1 *1/*5 -> Decreased_Function -> Simvastatin alternative_therapy
                    {"sample_id": "PT01", "rsid": "rs116855232", "genotype": "C/T"}, # NUDT15 *1/*3 -> IM -> Thiopurine dose_reduction_required
                    {"sample_id": "PT01", "rsid": "rs887829", "genotype": "T/T"},    # UGT1A1 *80/*80 -> PM -> Irinotecan dose_reduction_required
                ]
            ).to_csv(pgx_path, sep="\t", index=False)

            res = _run_json(
                "germ_pgx_star_allele.py",
                ["--genotypes", str(pgx_path), "--output_dir", str(tdir / "pgx_out")],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "pgx-star-allele-result.json")
            self.assertEqual(res["input_mode"], "rsid_table")
            self.assertGreater(res["n_actionable_alerts"], 0)

            pdf = pd.read_csv(res["outputs"]["cpic_diplotype_phenotypes"], sep="\t").set_index("gene")
            self.assertEqual(pdf.loc["CYP2C19", "diplotype"], "*2/*2")
            self.assertEqual(pdf.loc["CYP2C19", "phenotype"], "Poor_Metabolizer")
            self.assertEqual(pdf.loc["DPYD", "diplotype"], "*1/*2A")
            self.assertEqual(float(pdf.loc["DPYD", "activity_score"]), 1.0)
            self.assertEqual(pdf.loc["DPYD", "phenotype"], "Intermediate_Metabolizer")

            rdf = pd.read_csv(res["outputs"]["cpic_drug_recommendations"], sep="\t")
            clop = rdf[(rdf["gene"] == "CYP2C19") & (rdf["drug"] == "Clopidogrel")].iloc[0]
            self.assertEqual(clop["action_tier"], "alternative_therapy")
            fu5 = rdf[rdf["gene"] == "DPYD"].iloc[0]
            self.assertEqual(fu5["action_tier"], "dose_reduction_required")
            cod = rdf[(rdf["gene"] == "CYP2D6") & (rdf["drug"] == "Codeine")].iloc[0]
            self.assertEqual(cod["action_tier"], "contraindicated_or_avoid")

    def test_5_cnv_sv_integrate_snv_del_compound(self) -> None:
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            sv_path = tdir / "germline_svs.tsv"
            pd.DataFrame(
                [
                    # 1) Exon deletion in AR gene PAH (heterozygous CN=1) -> pairs with SNV!
                    {"sample": "PT_PKU", "chrom": "chr12", "start": 102850000, "end": 102865000, "sv_type": "DEL", "copy_number": 1, "genes": "PAH"},
                    # 2) PMP22 triplosensitive duplication -> Pathogenic
                    {"sample": "PT_PKU", "chrom": "chr17", "start": 15225000, "end": 15270000, "sv_type": "DUP", "copy_number": 3, "genes": "PMP22"},
                    # 3) DMD haploinsufficient deletion -> Pathogenic
                    {"sample": "PT_PKU", "chrom": "chrX", "start": 32000000, "end": 32150000, "sv_type": "DEL", "copy_number": 0, "genes": "DMD"},
                ]
            ).to_csv(sv_path, sep="\t", index=False)

            snv_path = tdir / "snv_candidates.tsv"
            pd.DataFrame(
                [
                    {
                        "sample": "PT_PKU",
                        "variant_id": "chr12:102840000:C>T",
                        "gene": "PAH",
                        "consequence": "missense_variant",
                        "proband_gt": "0/1",
                        "gnomad_af": 0.0005,
                        "acmg_class": "Likely_pathogenic",
                    }
                ]
            ).to_csv(snv_path, sep="\t", index=False)

            res = _run_json(
                "germ_cnv_sv_integrate.py",
                [
                    "--sv_calls", str(sv_path),
                    "--snv_candidates", str(snv_path),
                    "--output_dir", str(tdir / "cnv_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res, "cnv-sv-integrate-result.json")
            self.assertEqual(res["n_sv_calls"], 3)
            self.assertEqual(res["n_snv_cnv_compound_hits"], 1)
            self.assertEqual(res["compound_hits"][0]["gene"], "PAH")
            self.assertEqual(res["compound_hits"][0]["compound_mechanism"], "trans_snv_plus_del_compound_het")

            sv_df = pd.read_csv(res["outputs"]["classified_germline_cnv_sv"], sep="\t").set_index("genes")
            self.assertEqual(sv_df.loc["PMP22", "acmg_cnv_class"], "Pathogenic")
            self.assertEqual(sv_df.loc["DMD", "acmg_cnv_class"], "Pathogenic")
            self.assertEqual(sv_df.loc["PAH", "acmg_cnv_class"], "Likely_pathogenic")


if __name__ == "__main__":
    unittest.main()
