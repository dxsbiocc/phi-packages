# proteomics

Downstream quantitative MS-based proteomics and PTM/phosphoproteomics plugin
(`toolPrefix: prot`, private pure-conda environment `plugin:proteomics`).

## Scope & division of work

- **Upstream (out of scope here)**: raw `.raw` / `.d` / `.mzML` spectrum
  searching, FDR scoring, and quantification (DIA-NN, MaxQuant, FragPipe,
  Spectronaut, or `nf-core/quantms` wrappers).
- **This plugin (`prot_*`)**:
  1. `prot_rollup` — aggregate peptide/precursor intensity tables (DIA-NN
     `pr_matrix`, MaxQuant `peptides.txt`, FragPipe `combined_peptide`, or
     QUANTIX/generic peptide tables) to log2 protein abundances (`rrollup`,
     `top3`, `median`, `sum`) with `unique_only` or `razor` shared-peptide
     assignment.
  2. `prot_inspect` — auto-detect DIA-NN / MaxQuant / FragPipe / generic protein
     tables, contaminant/reverse rows, scale (`raw_intensity` vs `log2`),
     per-sample missingness, sample alignment, replicate counts, and
     batch-condition collinearity.
  2. `prot_qc_impute` — contaminant/decoy filtering, log2 transformation,
     group-aware valid-value filtering, median/quantile normalization,
     Perseus-style left-censored MinProb (MNAR) or kNN (MAR) imputation, and
     PCA / CV / missingness QC plots.
  3. `prot_de` — empirical-Bayes moderated t-test (limma-style variance
     moderation), Welch t-test, or Mann-Whitney U with BH FDR, volcano and MA
     plots.
  4. `prot_ptm` — PTM/phosphosite differential analysis with optional
     total-protein abundance normalization (`site_log2 - protein_log2`) and
     Kinase-Substrate Enrichment Analysis (KSEA).
  5. `prot_enrichment` — offline hypergeometric ORA (strictly against the
     detected-proteome background universe) and deterministic pre-ranked GSEA
     against GMT gene sets.
  6. `prot_network` — co-abundance correlation network, connected module
     detection, weighted hub scoring, and clustered correlation heatmap.
