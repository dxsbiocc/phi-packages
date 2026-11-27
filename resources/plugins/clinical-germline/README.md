# clinical-germline · 0.1.0

Downstream rare-disease germline variant interpretation and clinical pharmacogenomics plugin (`toolPrefix: germ`, agent `ClinicalGermline`) running in a locked MIT/BSD `conda-forge` environment (`plugin:germ`: Python 3.13 + `pandas`/`numpy`/`scipy`/`scikit-learn`/`matplotlib-base`).

All code in `plugin:germ` is 100% original MIT-licensed Python with zero non-commercial `InterVar` or `ANNOVAR` source code, implementing published ACMG/AMP 2015 (Richards et al.), Tavtigian et al. (2018/2020) Bayesian point-scale rules, ClinGen/Pejaver (2022) calibrated computational thresholds, ClinGen/ACMG (2020) Riggs dosage-sensitivity scoring, and CPIC diplotype-to-phenotype guidelines.

## Division of Work: Upstream `nf-core` Wrappers & Connectors vs. `germ_*` Tools

| Stage / Task | Upstream Wrapper / Connector | Deterministic Downstream Tool (`plugin:germ`) |
| --- | --- | --- |
| End-to-end germline WGS/WES alignment & calling | `nf-core/sarek` | — |
| Deep learning germline SNV/indel calling | `nf-core/modules/deepvariant` | — |
| Long-read / hybrid germline variant calling | `nf-core/modules/clair3` | — |
| Germline structural variant (SV) calling | `nf-core/modules/manta` | — |
| Functional transcript annotation (CSQ / ANN) | `nf-core/modules/ensemblvep`, `nf-core/modules/snpeff` | — |
| Interactive HPO ontology lookup & semantic expansion | `ols4` MCP connector (EMBL-EBI OLS4) | — |
| Trio / singleton Mendelian inheritance & phasing filter | — | `germ_vcf_filter` |
| ACMG/AMP 2015 + Tavtigian Bayesian & SF v3.2 classification | — | `germ_acmg_prioritize` |
| IC-weighted HPO phenotype similarity & diagnostic ranking | — | `germ_phenotype_rank` |
| CPIC star-allele calling, Activity Scores & Level A dosing | — | `germ_pgx_star_allele` |
| ClinGen CNV/SV dosage scoring & SNV+DEL compound hits | — | `germ_cnv_sv_integrate` |

## Tools (`germ_*`)

- `germ_vcf_filter`: Parses single-sample or trio VCF/VCF.gz (`INFO/CSQ`, `INFO/ANN`, `FORMAT/GT:AD:DP:GQ`) or annotated TSV/CSV tables; applies depth/GQ and coding/splice gates; evaluates `de_novo`, `autosomal_recessive_hom`, `compound_heterozygous` (*in trans* vs `unphased_candidate_comp_het` in singletons), `autosomal_dominant`, and `x_linked` fits; outputs `filtered_germline_candidates.tsv`, `inheritance_gene_summary.tsv`, and `germline_inheritance_funnel.png`.
- `germ_acmg_prioritize`: Evaluates ACMG/AMP 2015 criteria (`PVS1`, `PS1`/`PS5_ClinVar`, `PS2`, `PM2_Supporting`, `PM3`, `PM4`, Pejaver-calibrated `PP3`, `BA1`, `BS1`, `BP4`, `BP7`) under the Tavtigian Bayesian point scale ($P_0 = 0.10$, $O_{\text{PVS}} = 350$), flags ACMG SF v3.2 actionable secondary findings, and outputs `acmg_classified_variants.tsv` and `acmg_evidence_summary.png`.
- `germ_phenotype_rank`: Computes Information-Content-weighted phenotypic similarity (`phenotype_sim` in `0..1`) between patient HPO terms/keywords (`HP:XXXXXXX` or English/Chinese symptoms) and a built-in ~60-gene rare-disease HPO map (or custom `genes_to_phenotype.txt`), integrating ACMG points and inheritance fit into an Exomiser-inspired `combined_diagnostic_score` (`0..100`).
- `germ_pgx_star_allele`: Calls/resolves star-allele diplotypes (`CYP2C19`, `CYP2D6`, `CYP2C9`, `VKORC1`, `SLCO1B1`, `DPYD`, `TPMT`, `NUDT15`, `UGT1A1`), CPIC Activity Scores, metabolizer/transporter phenotypes, and CPIC Level A clinical dosing recommendations (`contraindicated_or_avoid`, `dose_reduction_required`, `alternative_therapy`, `standard_dosing`).
- `germ_cnv_sv_integrate`: Scores germline CNV/SV calls (`DEL`, `DUP`, `INV`, `BND`) using ClinGen/ACMG 2020 Riggs dosage-sensitivity criteria and cross-references rare coding SNVs to detect *trans*-acting **SNV + exon-DEL** compound heterozygotes in recessive disease genes.
