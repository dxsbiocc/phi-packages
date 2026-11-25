#!/usr/bin/env python3
"""End-to-end tests for pure-stdlib P2 upstream/downstream bridging tools:
- gen_vcf2maf (genomics 0.2.0): VEP CSQ / SnpEff ANN / ANNOVAR VCF -> GDC MAF + FASTA 3-mer Context
- rna_id_map (bulk-rna-seq 0.2.0): Ensembl version stripping, GTF biotype & symbol mapping, collision aggregation
- rna_geo_import (bulk-rna-seq 0.2.0): NCBI GEO GSE*_series_matrix.txt.gz + GPL probe collapse
"""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VCF2MAF = ROOT / "resources/plugins/genomics/skills/genomics-workflow/scripts/gen_vcf2maf.py"
RNA_ID_MAP = ROOT / "resources/plugins/bulk-rna-seq/skills/rnaseq-workflow/scripts/rna_id_map.py"
RNA_GEO_IMPORT = ROOT / "resources/plugins/bulk-rna-seq/skills/rnaseq-workflow/scripts/rna_geo_import.py"


def run_json(script: Path, args: list[str], cwd: Path) -> dict:
    proc = subprocess.run(
        [sys.executable, str(script), *args],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"{script.name} failed (rc={proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


class TestGenVcf2Maf(unittest.TestCase):
    def test_vep_snpeff_annovar_vcf_gz_with_fasta_context_and_sidecar(self):
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # Synthetic chr17 FASTA (30 bp: 1-based pos 10 has triplet A C G)
            fasta = tdir / "ref.fa"
            # pos 1..9 = TTTTTTTTT, pos 9..11 = A C G, pos 12..30 = T...
            seq = "TTTTTTTTACGTTTTTTTTTTTTTTTTTTT"
            fasta.write_text(f">chr17\n{seq}\n", encoding="utf-8")

            vcf_gz = tdir / "somatic.vcf.gz"
            vcf_lines = [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|HGVSc|HGVSp">',
                '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: \'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p\'">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR_01\tNORMAL_01",
                # 1) VEP multi-transcript SNP at pos 10 (ref C -> alt T), context at pos 9..11 = ACG; normal is 0/0
                "chr17\t10\t.\tC\tT\t.\tPASS\tCSQ=T|intron_variant|MODIFIER|TP53|ENSG1|Transcript|ENST_INTRON|protein_coding|c.100+5G>A|,T|missense_variant|MODERATE|TP53|ENSG1|Transcript|ENST_CANON|protein_coding|c.818G>A|p.Arg273His\tGT:AD:AF:DP\t0/1:60,40:0.40:100\t0/0:95,0:0.00:95",
                # 2) SnpEff frameshift deletion at pos 15 (ref TTG -> alt T => DEL of TG at 16..17)
                "chr17\t15\t.\tTTG\tT\t.\tPASS\tANN=T|frameshift_variant|HIGH|BRCA1|ENSG2|transcript|ENST_BRCA1|protein_coding|10/22|c.500_501delCA|p.Gln167fs\tGT:AD:DP\t0/1:50,30:80\t0/0:85,0:85",
                # 3) ANNOVAR stopgain SNP at pos 20 (ref T -> alt A)
                "chr17\t20\t.\tT\tA\t.\tPASS\tFunc.refGene=exonic;Gene.refGene=KRAS;ExonicFunc.refGene=stopgain;AAChange.refGene=KRAS:NM_004985:exon2:c.35T>A:p.L12X\tGT:AD:DP\t0/1:70,30:100\t0/0:90,0:90",
                # 4) LowQual record that must be dropped by default and kept with --include_filtered
                "chr17\t25\t.\tT\tG\t.\tLowQual\tCSQ=G|synonymous_variant|LOW|EGFR|ENSG3|Transcript|ENST_EGFR|protein_coding|c.90T>G|p.Leu30=\tGT:AD:DP\t0/1:88,12:100\t0/0:92,0:92",
            ]
            with gzip.open(vcf_gz, "wt", encoding="utf-8") as fh:
                fh.write("\n".join(vcf_lines) + "\n")

            out_maf = tdir / "cohort.maf"
            res = run_json(
                VCF2MAF,
                ["--vcf", str(vcf_gz), "--output_maf", str(out_maf), "--ref_fasta", str(fasta)],
                cwd=tdir,
            )
            self.assertEqual(res["n_vcf_records"], 4)
            self.assertEqual(res["n_filtered_out"], 1)
            self.assertEqual(res["n_maf_variants"], 3)
            self.assertEqual(res["tumor_samples"], ["TUMOR_01"])
            self.assertEqual(res["matched_normal"], "NORMAL_01")
            self.assertTrue(res["context_extracted"])
            self.assertEqual(sorted(res["annotators_detected"]), ["annovar", "snpeff", "vep"])

            with open(out_maf, "r", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(len(rows), 3)
            tp53 = rows[0]
            self.assertEqual(tp53["Hugo_Symbol"], "TP53")
            self.assertEqual(tp53["Variant_Classification"], "Missense_Mutation")
            self.assertEqual(tp53["Variant_Type"], "SNP")
            self.assertEqual(tp53["Context"], "ACG")
            self.assertAlmostEqual(float(tp53["t_vaf"]), 0.40, places=3)

            brca1 = rows[1]
            self.assertEqual(brca1["Hugo_Symbol"], "BRCA1")
            self.assertEqual(brca1["Variant_Classification"], "Frame_Shift_Del")
            self.assertEqual(brca1["Variant_Type"], "DEL")
            self.assertEqual(brca1["Start_Position"], "16")
            self.assertEqual(brca1["End_Position"], "17")
            self.assertEqual(brca1["Reference_Allele"], "TG")
            self.assertEqual(brca1["Tumor_Seq_Allele2"], "-")

            kras = rows[2]
            self.assertEqual(kras["Hugo_Symbol"], "KRAS")
            self.assertEqual(kras["Variant_Classification"], "Nonsense_Mutation")

            sidecar = Path(str(out_maf) + ".phi-artifact.json")
            self.assertTrue(sidecar.is_file())
            sc = json.loads(sidecar.read_text(encoding="utf-8"))
            self.assertEqual(sc["contractVersion"], "1.0.0")
            self.assertEqual(sc["kind"], "table")

            # Verify Phi boolean flag --include_filtered (no value argument) retains LowQual
            out_all = tdir / "cohort_all.maf"
            res_all = run_json(
                VCF2MAF,
                ["--vcf", str(vcf_gz), "--output_maf", str(out_all), "--include_filtered"],
                cwd=tdir,
            )
            self.assertEqual(res_all["n_maf_variants"], 4)


class TestRnaBridgingTools(unittest.TestCase):
    def test_rna_id_map_counts_and_de_deduplication(self):
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            gtf = tdir / "genes.gtf"
            gtf.write_text(
                "\n".join(
                    [
                        'chr17\tENSEMBL\tgene\t1\t1000\t.\t+\t.\tgene_id "ENSG00000141510.18"; gene_name "TP53"; gene_biotype "protein_coding";',
                        'chr17\tENSEMBL\tgene\t2000\t3000\t.\t+\t.\tgene_id "ENSG00000999999.2"; gene_name "TP53"; gene_biotype "protein_coding";',
                        'chr7\tENSEMBL\tgene\t4000\t5000\t.\t+\t.\tgene_id "ENSG00000146648.15"; gene_name "EGFR"; gene_biotype "protein_coding";',
                        'chr1\tENSEMBL\tgene\t6000\t7000\t.\t+\t.\tgene_id "ENSG00000222222.1"; gene_name "LINC001"; gene_biotype "lncRNA";',
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            counts_in = tdir / "raw_counts.tsv"
            counts_in.write_text(
                "gene_id\tS1\tS2\tS3\n"
                "ENSG00000141510.18\t100\t120\t140\n"
                "ENSG00000999999.2\t15\t20\t10\n"
                "ENSG00000146648.15\t50\t60\t70\n"
                "ENSG00000222222.1\t30\t40\t50\n",
                encoding="utf-8",
            )
            counts_out = tdir / "symbol_counts.tsv"
            res = run_json(
                RNA_ID_MAP,
                [
                    "--input",
                    str(counts_in),
                    "--output",
                    str(counts_out),
                    "--gtf",
                    str(gtf),
                    "--biotype",
                    "protein_coding",
                ],
                cwd=tdir,
            )
            self.assertEqual(res["table_kind"], "matrix")
            self.assertEqual(res["agg_method"], "sum")
            self.assertEqual(res["n_versions_stripped"], 4)
            self.assertEqual(res["n_biotype_filtered"], 1)
            self.assertEqual(res["n_collisions_merged"], 1)
            self.assertEqual(res["n_output_genes"], 2)
            self.assertEqual(res["scientificContext"][str(counts_out)]["valueSemantics"], "raw_counts")
            self.assertTrue(Path(res["collision_report"]).is_file())
            self.assertTrue(Path(str(counts_out) + ".phi-artifact.json").is_file())

            with open(counts_out, "r", encoding="utf-8") as fh:
                mapped = {r["gene_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(mapped["TP53"]["S1"], "115")
            self.assertEqual(mapped["TP53"]["S2"], "140")
            self.assertNotIn("LINC001", mapped)

            # DE table deduplication via min_pvalue
            de_in = tdir / "deseq2.tsv"
            de_in.write_text(
                "gene_id\tlog2FoldChange\tstat\tpvalue\tpadj\n"
                "ENSG00000141510.18\t2.5\t6.1\t1e-8\t1e-6\n"
                "ENSG00000999999.2\t0.2\t0.5\t0.62\t0.85\n"
                "ENSG00000146648.15\t-1.8\t-4.2\t2e-5\t4e-4\n",
                encoding="utf-8",
            )
            de_out = tdir / "deseq2_symbols.tsv"
            res_de = run_json(
                RNA_ID_MAP,
                ["--input", str(de_in), "--output", str(de_out), "--gtf", str(gtf)],
                cwd=tdir,
            )
            self.assertEqual(res_de["table_kind"], "de")
            self.assertEqual(res_de["agg_method"], "min_pvalue")
            self.assertEqual(res_de["scientificContext"][str(de_out)]["valueSemantics"], "de_stats")
            with open(de_out, "r", encoding="utf-8") as fh:
                de_rows = {r["gene_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(de_rows["TP53"]["log2FoldChange"], "2.5")

    def test_rna_geo_import_series_matrix_and_gpl_collapse(self):
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            series_gz = tdir / "GSE99999_series_matrix.txt.gz"
            lines = [
                '!Series_title\t"Synthetic Drug Response Study"',
                '!Series_geo_accession\t"GSE99999"',
                '!Series_platform_id\t"GPL570"',
                '!Sample_title\t"Ctrl_1"\t"Ctrl_2"\t"Treat_1"\t"Treat_2"',
                '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"',
                '!Sample_source_name_ch1\t"Lung"\t"Lung"\t"Lung"\t"Lung"',
                '!Sample_organism_ch1\t"Homo sapiens"\t"Homo sapiens"\t"Homo sapiens"\t"Homo sapiens"',
                '!Sample_characteristics_ch1\t"treatment: vehicle"\t"treatment: vehicle"\t"treatment: osimertinib"\t"treatment: osimertinib"',
                '!Sample_characteristics_ch1\t"age: 58"\t"age: 62"\t"age: 55"\t"age: 64"',
                "!series_matrix_table_begin",
                '"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"',
                '"1001_at"\t8.2\t8.4\t10.6\t10.8',
                '"2001_at"\t5.1\t5.0\t4.9\t5.2',
                '"2002_at"\t9.5\t9.7\t6.8\t6.5',
                '"3001_at"\t4.0\t4.1\t4.2\t4.0',
                "!series_matrix_table_end",
            ]
            with gzip.open(series_gz, "wt", encoding="utf-8") as fh:
                fh.write("\n".join(lines) + "\n")

            gpl = tdir / "GPL570_annot.tsv"
            gpl.write_text(
                "ID\tGene Symbol\n"
                "1001_at\tTP53 /// WRAP53\n"
                "2001_at\tEGFR\n"
                "2002_at\tEGFR\n"
                "3001_at\t---\n",
                encoding="utf-8",
            )

            outdir = tdir / "geo_out"
            res = run_json(
                RNA_GEO_IMPORT,
                [
                    "--series_matrix",
                    str(series_gz),
                    "--platform_table",
                    str(gpl),
                    "--output_dir",
                    str(outdir),
                ],
                cwd=tdir,
            )
            self.assertEqual(res["series_accession"], "GSE99999")
            self.assertEqual(res["n_samples"], 4)
            self.assertIn("treatment", res["metadata_columns"])
            self.assertIn("age", res["metadata_columns"])
            self.assertEqual(res["n_features_raw"], 4)
            self.assertEqual(res["n_features_output"], 2)
            self.assertEqual(res["inferred_scale"], "log_transformed")

            mat_path = Path(res["outputs"]["expression_matrix"])
            ss_path = Path(res["outputs"]["samplesheet"])
            self.assertTrue(mat_path.is_file())
            self.assertTrue(ss_path.is_file())
            self.assertTrue(Path(str(mat_path) + ".phi-artifact.json").is_file())
            self.assertTrue(Path(str(ss_path) + ".phi-artifact.json").is_file())

            with open(mat_path, "r", encoding="utf-8") as fh:
                mat_rows = {r["gene_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
            # 1001_at maps to primary symbol TP53; EGFR collapses to higher-mean probe 2002_at (9.5 vs 5.1)
            self.assertIn("TP53", mat_rows)
            self.assertEqual(mat_rows["EGFR"]["GSM1"], "9.5")


if __name__ == "__main__":
    unittest.main()
