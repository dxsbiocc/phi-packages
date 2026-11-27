#!/usr/bin/env python3
"""Biomni-Eval1 style cross-plugin end-to-end benchmark chains (4 biological workflows)."""
from __future__ import annotations

import importlib.util
import json
import math
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLUGINS_DIR = ROOT / "resources" / "plugins"

HAS_SCIPY_STACK = all(
    importlib.util.find_spec(pkg) is not None
    for pkg in ("numpy", "pandas", "scipy", "sklearn", "matplotlib")
)


def _run_json(script_path: Path, args: list[str], cwd: Path) -> dict:
    proc = subprocess.run(
        [sys.executable, "-B", str(script_path), *args],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"{script_path.name} failed (rc={proc.returncode}):\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


@unittest.skipUnless(HAS_SCIPY_STACK, "numpy/pandas/scipy/sklearn/matplotlib required for live E2E benchmarks")
class TestBiomniEvalBenchmarks(unittest.TestCase):
    def test_chain_1_gwas_finemap_to_mendelian_randomization(self) -> None:
        """Biomni-Eval Chain 1: statistical-genetics (gwas_sumstats_harmonize + gwas_finemap) -> mendelian-randomization (mr_harmonize + mr_estimate)."""
        import numpy as np
        import pandas as pd
        from scipy import stats as sp_stats

        rng = random.Random(2026)
        gwas_scripts = PLUGINS_DIR / "statistical-genetics" / "skills" / "gwas-workflow" / "scripts"
        mr_scripts = PLUGINS_DIR / "mendelian-randomization" / "skills" / "mr-workflow" / "scripts"

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            exp_rows = []
            out_rows = []
            true_causal = 0.45

            # 18 strong non-palindromic instruments + 1 flipped + 1 ambiguous palindromic + 2 null/weak
            for i in range(22):
                snp = f"rs{10000 + i}"
                chrom = 1 + (i // 8)
                pos = 1_000_000 + (i % 8) * 700_000
                if i < 18:
                    ea_e, nea_e = "A", "G"
                    eaf = round(0.18 + 0.015 * i, 3)
                    bx = 0.08 + 0.008 * i
                    sex = 0.011
                    ea_o, nea_o, eaf_o, sign_o = ea_e, nea_e, eaf, 1.0
                    if i == 5:
                        ea_o, nea_o, eaf_o, sign_o = nea_e, ea_e, round(1.0 - eaf, 3), -1.0
                elif i == 18:
                    ea_e, nea_e, eaf = "A", "T", 0.22
                    bx, sex = 0.12, 0.012
                    ea_o, nea_o, eaf_o, sign_o = "T", "A", 0.78, -1.0
                elif i == 19:
                    ea_e, nea_e, eaf = "C", "G", 0.49
                    bx, sex = 0.14, 0.012
                    ea_o, nea_o, eaf_o, sign_o = "C", "G", 0.49, 1.0
                else:
                    ea_e, nea_e, eaf = "C", "T", 0.30
                    bx, sex = 0.01, 0.02
                    ea_o, nea_o, eaf_o, sign_o = "C", "T", 0.30, 1.0

                px = float(2.0 * sp_stats.norm.sf(abs(bx / sex)))
                by = sign_o * (true_causal * bx + rng.gauss(0.0, 0.004))
                sey = 0.015
                py = float(2.0 * sp_stats.norm.sf(abs(by / sey)))
                exp_rows.append(
                    {"ID": snp, "CHROM": chrom, "GENPOS": pos, "ALLELE0": nea_e, "ALLELE1": ea_e, "A1FREQ": eaf, "INFO": 0.96, "BETA": bx, "SE": sex, "P": px, "N": 50000}
                )
                out_rows.append(
                    {"SNP": snp, "CHR": chrom, "POS": pos, "EA": ea_o, "NEA": nea_o, "EAF": eaf_o, "BETA": by, "SE": sey, "P": py, "N": 60000}
                )

            # Inject 2 low-QC rows that gwas_sumstats_harmonize must drop before MR
            exp_rows.append({"ID": "rs_low_maf", "CHROM": 1, "GENPOS": 9_500_000, "ALLELE0": "A", "ALLELE1": "G", "A1FREQ": 0.002, "INFO": 0.96, "BETA": 0.15, "SE": 0.02, "P": 1e-12, "N": 50000})
            exp_rows.append({"ID": "rs_low_info", "CHROM": 1, "GENPOS": 9_900_000, "ALLELE0": "C", "ALLELE1": "T", "A1FREQ": 0.25, "INFO": 0.40, "BETA": 0.15, "SE": 0.02, "P": 1e-12, "N": 50000})

            exp_raw = tdir / "exposure_regenie.tsv"
            out_raw = tdir / "outcome_gwas.tsv"
            pd.DataFrame(exp_rows).to_csv(exp_raw, sep="\t", index=False)
            pd.DataFrame(out_rows).to_csv(out_raw, sep="\t", index=False)

            harm_gwas = _run_json(
                gwas_scripts / "gwas_sumstats_harmonize.py",
                [
                    "--sumstats", str(exp_raw),
                    "--min_maf", "0.01",
                    "--min_info", "0.80",
                    "--lead_window_kb", "500",
                    "--output_dir", str(tdir / "step1_gwas_harmonize"),
                ],
                cwd=tdir,
            )
            self.assertEqual(harm_gwas["n_variants_harmonized"], 22)
            self.assertEqual(harm_gwas["n_dropped_qc"], 2)
            harm_exp_tsv = Path(harm_gwas["outputs"]["harmonized_sumstats"])
            self.assertTrue(harm_exp_tsv.is_file())

            # Step 2: Regional SuSiE-RSS fine-mapping of a GWAS lead locus
            m = 60
            n_fm = 25000.0
            r_mat = np.eye(m, dtype=float)
            for b0 in (0, 30):
                for i in range(b0, b0 + 30):
                    for j in range(b0, b0 + 30):
                        r_mat[i, j] = 0.85 ** abs(i - j)
            b_true = np.zeros(m, dtype=float)
            b_true[10] = 0.09
            b_true[45] = -0.08
            z_fm = math.sqrt(n_fm) * (r_mat @ b_true)
            p_fm = np.maximum(2.0 * sp_stats.norm.sf(np.abs(z_fm)), 1e-300)
            locus_path = tdir / "locus.tsv"
            ld_path = tdir / "locus_ld.txt"
            pd.DataFrame(
                {"SNP": [f"rs_fm_{i}" for i in range(m)], "CHR": ["6"] * m, "POS": [30_000_000 + i * 5_000 for i in range(m)], "Z": z_fm, "P": p_fm, "N": [n_fm] * m}
            ).to_csv(locus_path, sep="\t", index=False)
            np.savetxt(ld_path, r_mat, fmt="%.6f")

            fm_res = _run_json(
                gwas_scripts / "gwas_finemap.py",
                [
                    "--locus_sumstats", str(locus_path),
                    "--ld_matrix", str(ld_path),
                    "--max_causal", "5",
                    "--coverage", "0.95",
                    "--output_dir", str(tdir / "step2_finemap"),
                ],
                cwd=tdir,
            )
            self.assertEqual(fm_res["method"], "SuSiE-RSS")
            self.assertEqual(fm_res["n_credible_sets"], 2)

            # Step 3 & 4: Feed statistical-genetics harmonized summary stats directly into MR
            mr_harm = _run_json(
                mr_scripts / "mr_harmonize.py",
                [
                    "--exposure", str(harm_exp_tsv),
                    "--outcome", str(out_raw),
                    "--clump_kb", "500",
                    "--output_dir", str(tdir / "step3_mr_harm"),
                ],
                cwd=tdir,
            )
            self.assertGreaterEqual(mr_harm["n_harmonized_retained"], 18)
            self.assertGreaterEqual(mr_harm["n_palindromic_dropped"], 1)

            mr_est = _run_json(
                mr_scripts / "mr_estimate.py",
                [
                    "--harmonized", mr_harm["outputs"]["harmonized_table"],
                    "--output_dir", str(tdir / "step4_mr_est"),
                ],
                cwd=tdir,
            )
            ivw = mr_est["methods"][0]
            self.assertAlmostEqual(ivw["beta"], true_causal, delta=0.08)
            self.assertLess(ivw["pvalue"], 1e-6)
            self.assertTrue(mr_est["steiger_directionality"]["correct_causal_direction"])

    def test_chain_2_clinical_germline_trio_to_pheno_and_cnv_compound(self) -> None:
        """Biomni-Eval Chain 2: clinical-germline end-to-end diagnostic cascade (vcf_filter -> acmg_prioritize -> phenotype_rank -> cnv_sv_integrate)."""
        import pandas as pd

        germ_scripts = PLUGINS_DIR / "clinical-germline" / "skills" / "germline-workflow" / "scripts"
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            vcf_path = tdir / "proband_trio.vcf"
            vcf_lines = [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|HGVSc|HGVSp|gnomAD_AF|CLIN_SIG|REVEL|CADD_PHRED|SpliceAI_max|pLI|LOEUF">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND\tMOTHER\tFATHER",
                "chr2\t166053000\t.\tG\tA\t999\tPASS\tCSQ=A|stop_gained|HIGH|SCN1A|c.1234C>T|p.Arg412*|0.0|pathogenic||38.0|0.02|1.0|0.12\tGT:AD:DP:GQ\t0/1:24,22:46:99\t0/0:40,0:40:99\t0/0:44,0:44:99",
                "chr17\t43091000\t.\tCT\tC\t950\tPASS\tCSQ=C|frameshift_variant|HIGH|BRCA1|c.5266dup|p.Gln1756Profs|0.00002|pathogenic||35.0|0.0|0.99|0.18\tGT:AD:DP:GQ\t0/1:23,21:44:99\t0/1:25,20:45:99\t0/0:41,0:41:99",
                "chr12\t102840000\t.\tC\tT\t860\tPASS\tCSQ=T|missense_variant|MODERATE|PAH|c.1222C>T|p.Arg408Trp|0.00005|pathogenic|0.92|28.0|0.0|0.02|0.55\tGT:AD:DP:GQ\t0/1:22,20:42:99\t0/0:38,0:38:99\t0/1:19,21:40:99",
                "chr1\t100000\t.\tA\tG\t50\tPASS\tCSQ=G|synonymous_variant|LOW|NOC2L|c.30A>G|p.Leu10=|0.00001||0.10|5.0|0.0|0.1|0.8\tGT:AD:DP:GQ\t0/1:3,2:5:12\t0/0:20,0:20:80\t0/0:20,0:20:80",
            ]
            vcf_path.write_text("\n".join(vcf_lines) + "\n", encoding="utf-8")

            step1 = _run_json(
                germ_scripts / "germ_vcf_filter.py",
                [
                    "--variants", str(vcf_path),
                    "--proband_id", "PROBAND",
                    "--mother_id", "MOTHER",
                    "--father_id", "FATHER",
                    "--output_dir", str(tdir / "step1_filter"),
                ],
                cwd=tdir,
            )
            self.assertTrue(step1["trio_mode"])

            step2 = _run_json(
                germ_scripts / "germ_acmg_prioritize.py",
                [
                    "--variants", step1["outputs"]["filtered_germline_candidates"],
                    "--output_dir", str(tdir / "step2_acmg"),
                ],
                cwd=tdir,
            )
            self.assertGreaterEqual(step2["n_pathogenic_or_likely"], 2)
            self.assertEqual(step2["n_acmg_sf_reportable"], 1)

            step3 = _run_json(
                germ_scripts / "germ_phenotype_rank.py",
                [
                    "--variants", step2["outputs"]["acmg_classified_variants"],
                    "--hpo_terms", "HP:0001250", "HP:0002373", "HP:0001263",
                    "--output_dir", str(tdir / "step3_pheno"),
                ],
                cwd=tdir,
            )
            rank_df = pd.read_csv(step3["outputs"]["diagnostic_gene_ranking"], sep="\t")
            self.assertEqual(rank_df.iloc[0]["gene"], "SCN1A")
            self.assertGreater(float(rank_df.iloc[0]["combined_diagnostic_score"]), 85.0)

            sv_path = tdir / "germline_svs.tsv"
            pd.DataFrame(
                [
                    {"sample": "PROBAND", "chrom": "chr12", "start": 102850000, "end": 102865000, "sv_type": "DEL", "copy_number": 1, "genes": "PAH"},
                    {"sample": "PROBAND", "chrom": "chr17", "start": 15225000, "end": 15270000, "sv_type": "DUP", "copy_number": 3, "genes": "PMP22"},
                ]
            ).to_csv(sv_path, sep="\t", index=False)

            step4 = _run_json(
                germ_scripts / "germ_cnv_sv_integrate.py",
                [
                    "--sv_calls", str(sv_path),
                    "--snv_candidates", step2["outputs"]["acmg_classified_variants"],
                    "--output_dir", str(tdir / "step4_cnv_sv"),
                ],
                cwd=tdir,
            )
            self.assertEqual(step4["n_snv_cnv_compound_hits"], 1)
            self.assertEqual(step4["compound_hits"][0]["gene"], "PAH")
            self.assertEqual(step4["compound_hits"][0]["compound_mechanism"], "trans_snv_plus_del_compound_het")

    def test_chain_3_cancer_genomics_to_immunoinformatics(self) -> None:
        """Biomni-Eval Chain 3: genomics (gen_hrd_msi + gen_clonality) -> immunoinformatics (imm_hla_summary + imm_neoantigen_prioritize)."""
        import pandas as pd

        gen_scripts = PLUGINS_DIR / "genomics" / "skills" / "genomics-workflow" / "scripts"
        imm_scripts = PLUGINS_DIR / "immunoinformatics" / "skills" / "immuno-workflow" / "scripts"

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # Step 1: HRD/MSI genomic scar + microsatellite biomarker stratification
            seg_path = tdir / "ascat_segments.tsv"
            pd.DataFrame(
                [
                    {"sample": "TUMOR_A", "chrom": "chr1", "start": 1, "end": 22000000, "nMajor": 2, "nMinor": 0},
                    {"sample": "TUMOR_A", "chrom": "chr1", "start": 25000000, "end": 43000000, "nMajor": 2, "nMinor": 1},
                    {"sample": "TUMOR_A", "chrom": "chr1", "start": 43200001, "end": 60000000, "nMajor": 3, "nMinor": 1},
                    {"sample": "TUMOR_A", "chrom": "chr1", "start": 220000000, "end": 248500000, "nMajor": 3, "nMinor": 1},
                    {"sample": "TUMOR_A", "chrom": "chr3", "start": 120000000, "end": 145000000, "nMajor": 2, "nMinor": 0},
                    {"sample": "TUMOR_B", "chrom": "chr1", "start": 1, "end": 248000000, "nMajor": 1, "nMinor": 1},
                ]
            ).to_csv(seg_path, sep="\t", index=False)

            msi_path = tdir / "msisensor.tsv"
            pd.DataFrame(
                [
                    {"sample": "TUMOR_A", "Total_Number_of_Sites": 1000, "Number_of_Somatic_Sites": 8, "%": 0.8},
                    {"sample": "TUMOR_B", "Total_Number_of_Sites": 1000, "Number_of_Somatic_Sites": 210, "%": 21.0},
                ]
            ).to_csv(msi_path, sep="\t", index=False)

            hrd_res = _run_json(
                gen_scripts / "gen_hrd_msi.py",
                [
                    "--segments", str(seg_path),
                    "--msi_table", str(msi_path),
                    "--genome", "hg38",
                    "--hrd_cutoff", "5",
                    "--msi_pct_cutoff", "3.5",
                    "--output_dir", str(tdir / "step1_hrd_msi"),
                ],
                cwd=tdir,
            )
            by_s = {r["sample"]: r for r in hrd_res["samples"]}
            self.assertEqual(by_s["TUMOR_A"]["biomarker_stratum"], "HRD_pos_MSS")
            self.assertEqual(by_s["TUMOR_B"]["biomarker_stratum"], "MSI_H")

            # Step 2: Clonality / CCF inference for TUMOR_B (MSI-H checkpoint candidate)
            pur_path = tdir / "purity.tsv"
            pd.DataFrame([{"sample": "TUMOR_B", "purity": 0.60, "ploidy": 2.0}]).to_csv(pur_path, sep="\t", index=False)
            maf_path = tdir / "tumor_b.maf"
            pd.DataFrame(
                [
                    {"Hugo_Symbol": "KRAS", "Tumor_Sample_Barcode": "TUMOR_B", "Chromosome": "12", "Start_Position": 25398284, "Variant_Classification": "Missense_Mutation", "t_alt_count": 30, "t_ref_count": 70},
                    {"Hugo_Symbol": "TP53", "Tumor_Sample_Barcode": "TUMOR_B", "Chromosome": "17", "Start_Position": 7577120, "Variant_Classification": "Missense_Mutation", "t_alt_count": 29, "t_ref_count": 71},
                    {"Hugo_Symbol": "SUBCLONE_G", "Tumor_Sample_Barcode": "TUMOR_B", "Chromosome": "1", "Start_Position": 15000000, "Variant_Classification": "Missense_Mutation", "t_alt_count": 6, "t_ref_count": 94},
                ]
            ).to_csv(maf_path, sep="\t", index=False)

            clon_res = _run_json(
                gen_scripts / "gen_clonality.py",
                [
                    "--maf", str(maf_path),
                    "--purity_table", str(pur_path),
                    "--output_dir", str(tdir / "step2_clonality"),
                ],
                cwd=tdir,
            )
            ccf_df = pd.read_csv(clon_res["variants_table"], sep="\t").set_index("Hugo_Symbol")
            self.assertEqual(ccf_df.loc["KRAS", "clonality"], "clonal")
            self.assertEqual(ccf_df.loc["SUBCLONE_G", "clonality"], "subclonal")

            # Step 3: Consensus HLA genotyping for TUMOR_B
            opt_path = tdir / "TUMOR_B_result.tsv"
            opt_path.write_text(
                "\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n"
                "0\tA*02:01\tA*01:01\tB*07:02\tB*44:02\tC*04:01\tC*07:02\t5000\t4900\n",
                encoding="utf-8",
            )
            hla_res = _run_json(
                imm_scripts / "imm_hla_summary.py",
                ["--hla_inputs", str(opt_path), "--output_dir", str(tdir / "step3_hla")],
                cwd=tdir,
            )
            hla_row = pd.read_csv(hla_res["outputs"]["hla_consensus_genotypes"], sep="\t").set_index("sample_id").loc["TUMOR_B"]
            patient_a1 = str(hla_row["HLA_A_1"])

            # Step 4: Prioritize somatic neoantigens using CCF directly from gen_clonality + HLA allele from imm_hla_summary
            neo_cand = tdir / "neoantigen_candidates.tsv"
            pd.DataFrame(
                [
                    {
                        "gene": "KRAS",
                        "hla_allele": patient_a1,
                        "mt_epitope": "KLVVVGAVGV",
                        "wt_epitope": "KLVVVGAGGV",
                        "ic50_mt": 42.0,
                        "ic50_wt": 3800.0,
                        "tpm": 32.0,
                        "ccf": float(ccf_df.loc["KRAS", "ccf"]),
                    },
                    {
                        "gene": "SUBCLONE_G",
                        "hla_allele": patient_a1,
                        "mt_epitope": "YLQPRTFLL",
                        "wt_epitope": "YSQPRTFLL",
                        "ic50_mt": 65.0,
                        "ic50_wt": 2000.0,
                        "tpm": 0.1,
                        "ccf": float(ccf_df.loc["SUBCLONE_G", "ccf"]),
                    },
                ]
            ).to_csv(neo_cand, sep="\t", index=False)

            neo_res = _run_json(
                imm_scripts / "imm_neoantigen_prioritize.py",
                [
                    "--candidates", str(neo_cand),
                    "--ic50_cutoff_nm", "500",
                    "--min_tpm", "1.0",
                    "--output_dir", str(tdir / "step4_neoantigen"),
                ],
                cwd=tdir,
            )
            self.assertEqual(neo_res["n_passed_all_gates"], 1)
            ndf = pd.read_csv(neo_res["outputs"]["prioritized_neoantigens"], sep="\t").set_index("gene")
            self.assertTrue(bool(ndf.loc["KRAS", "passes_all_gates"]))
            self.assertFalse(bool(ndf.loc["SUBCLONE_G", "passes_expression_gate"]))

    def test_chain_4_protein_design_to_cloning_and_wetlab_spr(self) -> None:
        """Biomni-Eval Chain 4: protein-design (pd_developability_audit) -> molecular-cloning (mc_codon_optimize + mc_design_primers) -> wetlab-assay (wa_binding_kinetics)."""
        import numpy as np
        import pandas as pd

        pd_scripts = PLUGINS_DIR / "protein-design" / "skills" / "design-evaluation" / "scripts"
        mc_scripts = PLUGINS_DIR / "molecular-cloning" / "skills" / "cloning-workflow" / "scripts"
        wa_scripts = PLUGINS_DIR / "wetlab-assay" / "skills" / "assay-workflow" / "scripts"

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # Step 1: Developability triage selects clean soluble binder
            clean_aa = "MKLTVGDDAAKLTVGDDAAKLTVGDDAAWYA"
            risky_aa = "MAHCEDKAKHNGDGKNASTAEKKDALAALVALVALV"
            fa_path = tdir / "binders.fasta"
            fa_path.write_text(f">clean_binder\n{clean_aa}\n>risky_binder\n{risky_aa}\n", encoding="utf-8")

            audit_res = _run_json(
                pd_scripts / "developability_audit.py",
                ["--sequences_input", str(fa_path), "--buffer_ph", "7.4", "--output_dir", str(tdir / "step1_audit")],
                cwd=tdir,
            )
            self.assertEqual(audit_res["tier_counts"]["pass"], 1)
            self.assertEqual(audit_res["tier_counts"]["high_risk"], 1)

            seq_by_id = {"clean_binder": clean_aa, "risky_binder": risky_aa}
            audit_df = pd.read_csv(audit_res["outputs"]["developability_audit_tsv"], sep="\t")
            winner_row = audit_df[audit_df["developability_tier"] == "pass"].iloc[0]
            winner_id = str(winner_row["design_id"])
            self.assertEqual(winner_id, "clean_binder")

            # Step 2: E. coli codon optimize winning AA sequence -> GenBank construct
            winner_fa = tdir / "winner.fasta"
            winner_fa.write_text(f">{winner_id}\n{seq_by_id[winner_id]}\n", encoding="utf-8")
            opt_res = _run_json(
                mc_scripts / "mc_codon_optimize.py",
                [
                    "--sequence_input", str(winner_fa),
                    "--host", "E_coli_K12",
                    "--output_dir", str(tdir / "step2_codon"),
                ],
                cwd=tdir,
            )
            self.assertTrue(opt_res["translation_identical"])
            self.assertEqual(opt_res["remaining_restriction_sites"], 0)
            gb_path = opt_res["outputs"]["optimized_construct_gb"]
            self.assertTrue(Path(gb_path).is_file())

            # Step 3: Design GoldenGate cloning primers directly from GenBank template
            primer_res = _run_json(
                mc_scripts / "mc_design_primers.py",
                [
                    "--template", gb_path,
                    "--mode", "GoldenGate",
                    "--target_start", "1",
                    "--target_end", str(len(clean_aa) * 3),
                    "--goldengate_enzyme", "BsaI",
                    "--overhang_5p", "AATG",
                    "--overhang_3p", "GCTT",
                    "--output_dir", str(tdir / "step3_primers"),
                ],
                cwd=tdir,
            )
            self.assertGreater(primer_res["n_candidates"], 0)
            self.assertIn("GGTCTC", primer_res["top_pair"]["fwd_full"])

            # Step 4: SPR binding kinetics characterization of purified binder (true KD = 10 nM)
            rng = np.random.default_rng(2026)
            kon, koff, rmax, td = 2.0e5, 2.0e-3, 120.0, 900.0
            times = np.linspace(0.0, 1200.0, 121)
            spr_rows = []
            for c_nm in (1.25, 2.5, 5.0, 10.0, 20.0, 40.0):
                c_m = c_nm * 1e-9
                k_obs = kon * c_m + koff
                r_eq = (kon * c_m * rmax) / k_obs
                r_td = r_eq * (1.0 - math.exp(-k_obs * td))
                for t in times:
                    r = r_eq * (1.0 - math.exp(-k_obs * t)) if t <= td else r_td * math.exp(-koff * (t - td))
                    spr_rows.append({"time": round(float(t), 2), "conc_nm": c_nm, "response": round(r + float(rng.normal(0.0, 0.35)), 4)})
            spr_path = tdir / "binder_spr.tsv"
            pd.DataFrame(spr_rows).to_csv(spr_path, sep="\t", index=False)

            spr_res = _run_json(
                wa_scripts / "wa_binding_kinetics.py",
                [
                    "--sensorgrams", str(spr_path),
                    "--dissociation_start", "900",
                    "--conc_unit", "nM",
                    "--output_dir", str(tdir / "step4_spr"),
                ],
                cwd=tdir,
            )
            self.assertAlmostEqual(spr_res["kd_nm"], 10.0, delta=0.8)
            self.assertAlmostEqual(spr_res["kd_eq_nm"], 10.0, delta=2.5)


if __name__ == "__main__":
    unittest.main()
