# Changelog — clinical-germline

## 0.1.0 (2026-11-26)

Initial release of the rare-disease germline variant interpretation and clinical pharmacogenomics plugin (`toolPrefix: germ`, agent `ClinicalGermline`, private locked environment `plugin:germ`).

- `germ_vcf_filter`: Trio/singleton Mendelian inheritance filtering (`de_novo`, `autosomal_recessive_hom`, phased *in trans* `compound_heterozygous` vs singleton `unphased_candidate_comp_het`, `autosomal_dominant`, `x_linked`) from VCF (`CSQ`/`ANN`) or annotated TSV/CSV.
- `germ_acmg_prioritize`: Transparent ACMG/AMP 2015 + ClinGen Pejaver (2022) + Tavtigian (2020) Bayesian point-scale pathogenicity classification, curated 80+ gene LoF/HI baseline, and ACMG SF v3.2 secondary-finding flagging.
- `germ_phenotype_rank`: Information-Content-weighted HPO phenotype matching (`HP:XXXXXXX` and English/Chinese symptom keywords) with Exomiser-inspired `combined_diagnostic_score` (`0..100`) ranking.
- `germ_pgx_star_allele`: CPIC star-allele diplotype resolution, Activity Scores, metabolizer/transporter phenotypes, and Level A drug dosing alerts across `CYP2C19`, `CYP2D6`, `CYP2C9`, `VKORC1`, `SLCO1B1`, `DPYD`, `TPMT`, `NUDT15`, and `UGT1A1`.
- `germ_cnv_sv_integrate`: ClinGen/ACMG 2020 Riggs CNV/SV dosage-sensitivity classification and *trans*-acting SNV + exon-DEL compound heterozygote discovery.
