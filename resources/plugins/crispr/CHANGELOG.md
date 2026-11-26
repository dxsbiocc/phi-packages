# Changelog — crispr

## 0.1.0 (2026-11-26)

Initial release: CRISPR gene-editing design, off-target scoring, pooled functional screening, and post-editing outcome quantification plugin (`toolPrefix: crispr`, agent `Crispr`, private locked environment `plugin:crispr`).

- `crispr_design_guides`: multi-nuclease (`SpCas9`, `SpCas9-NG`, `SaCas9`, `Cas12a`) guide design with sequence-context efficiency scoring, Bae et al. (2014) MMEJ out-of-frame KO scoring, `CBE`/`ABE` base-editor window & bystander/stop-codon auditing, and Golden Gate cloning oligo synthesis.
- `crispr_offtarget_scan`: PAM-anchored mismatch off-target search (`0..4` mismatches, optional `NAG`/`NGA` non-canonical PAMs) with Hsu-Zhang **MIT** and Doench 2016 **CFD** specificity scoring.
- `crispr_library_check`: pooled sgRNA library pre-synthesis and pre-count QC (duplicate IDs/sequences, IUPAC checks, guides/gene, control representation, `TTTT` Pol III motifs, restriction enzyme collisions).
- `crispr_screen_qc`: pooled CRISPR screen QC (depth, zero-dropout %, Lorenz Gini index, DESeq median-of-ratios normalization, pairwise correlation & PCA, essential vs non-essential **NNMD** and ROC-AUC).
- `crispr_screen_hits`: $\alpha$-RRA gene-level positive/negative selection hit calling with moderated sgRNA variance, seeded empirical permutation calibration, Stouffer $Z$, and BH FDR.
- `crispr_editing_quant`: post-editing NHEJ frameshift vs in-frame KO, HDR / imperfect HDR, and `CBE`/`ABE` base-editor target vs bystander quantification from amplicon allele tables/FASTAs or `nf-core/modules/crispresso2` allele frequency tables.
