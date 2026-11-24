# Proteomics Downstream Analysis Playbook

## 1. Format & Scale Inspection (`prot_inspect`)
- **MaxQuant (`proteinGroups.txt`)**: look for `LFQ intensity ...` or `Intensity ...` columns and filter `Potential contaminant == "+"`, `Reverse == "+"`, `Only identified by site == "+"`.
- **DIA-NN (`report.pg_matrix.tsv`)**: protein group IDs in `Protein.Group` / `Genes`, sample columns are raw linear intensities (often > 1e4).
- **FragPipe (`combined_protein.tsv`)**: look for `... MaxLFQ Intensity` or `... Intensity` columns; `Protein` / `Gene` holds identifiers.
- **Scale check**: values with 95th percentile > 50 are almost certainly linear intensities and need `log2` transformation before normalization or differential testing.

## 2. Filtering, Normalization & Imputation (`prot_qc_impute`)
- **Group-aware valid-value filtering**: in MS proteomics, a protein genuinely off in control and on in treatment will have 50% overall missingness. Always prefer `filter_mode="in_at_least_one_group"` with `min_valid_frac=0.5`–`0.7` when sample metadata is available.
- **Normalization**: `median` centering (subtracting each sample's median log2 intensity and adding the global median) is the default standard for LFQ/DIA/TMT; use `quantile` only when intensity distributions differ strongly in shape across runs.
- **Imputation**:
  - `minprob` (Perseus-style left-censored Gaussian, `downshift=1.8`, `width=0.3`): best for LFQ/DIA where missing values are predominantly **MNAR** (below limit of detection).
  - `knn` (`knn_k=5`): best for TMT or shallow **MAR** missingness where dropouts are intensity-independent.
  - `none`: use when missingness is < 2% or when you prefer complete-case testing.

## 3. Differential Abundance (`prot_de`)
- **Empirical-Bayes moderated t-test (`moderated_t`, default)**: borrows variance information across all quantified proteins (Smyth 2004 limma framework), dramatically stabilizing small-n proteomics comparisons (`n = 3–6` per group).
- **Welch t-test (`welch`)** or **Mann-Whitney U (`mannwhitney`)**: available for sensitivity checks or large-n non-normal cohorts.
- Always inspect volcano and MA plots for global intensity-dependent bias before interpreting hits.

## 4. PTM / Phosphoproteomics (`prot_ptm`)
- A 2-fold increase in a phosphosite on a protein whose total abundance also doubled is **not** a stoichiometry change. Whenever a matched total-protein matrix exists, pass `--protein_matrix` so `prot_ptm` computes `site_log2 - protein_log2` before differential testing.
- **KSEA**: aggregates site-level log2FCs for each kinase's known substrates into a z-score `(mean_sub - mean_all) * sqrt(m) / sd_all` with two-sided p-value and BH FDR. Treat KSEA scores with `m < 3` substrates as unreliable.

## 5. Pathway Enrichment (`prot_enrichment`) & Co-Abundance (`prot_network`)
- **ORA background red line**: MS proteomics detects ~3,000–9,000 proteins, heavily biased toward abundant/soluble proteins. Testing against a 20,000-gene genome inflates false positives (ribosome, metabolism, cytoskeleton). Always pass the full `prot_de` table so the detected proteome serves as the universe.
- **GSEA**: threshold-free complement to ORA; ranks all quantified proteins by `moderated_t` (or signed `-log10(p)`).
- **Co-abundance network**: Spearman correlation across top variable/regulated proteins identifies co-regulated complexes and hub proteins (`strength = sum(|r|)` over significant edges).
