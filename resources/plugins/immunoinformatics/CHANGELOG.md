# Changelog — immunoinformatics

## 0.1.0 (2026-11-26)

Initial release of the downstream immunoinformatics plugin (`toolPrefix: imm`, agent `Immunoinformatics`, private locked environment `plugin:imm`).

- `imm_repertoire_inspect`: AIRR, 10x Cell Ranger V(D)J, TRUST4, and MiXCR harmonization, productivity/frame/stop QC, Kyte-Doolittle hydropathy, and CDR3 spectratype figure.
- `imm_diversity`: Shannon entropy, Pielou clonality, Gini-Simpson, Inverse Simpson, D50, clonal homeostasis tiers, V-J pairing matrix, and optional Mann-Whitney / Kruskal-Wallis + BH FDR comparisons.
- `imm_clonotype_tracking`: Pairwise Morisita-Horn and Jaccard overlap matrices, public clonotypes, enriched 4-mer CDR3 motifs, and Fisher's exact + BH FDR longitudinal clone expansion/contraction calling.
- `imm_hla_summary`: Multi-caller 2-field (`4-digit`) HLA consensus across Class I/II loci, Class I homozygosity immune-escape flagging, and approximate pocket-Grantham HLA Evolutionary Divergence (HED).
- `imm_neoantigen_prioritize`: MHC IC50 ingestion or anchor-motif PSSM estimation, Differential Agretopicity Index (DAI), anchor vs TCR-contact classification, foreignness/hydropathy, TPM/VAF/CCF gates, and composite immunogenicity scoring (`0..100`).
