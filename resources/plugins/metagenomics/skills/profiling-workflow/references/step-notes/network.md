# `mg_network` — Spearman co-occurrence network analysis

## Purpose

Constructs a taxon or functional-feature co-occurrence network from a
`features × samples` abundance matrix using pairwise Spearman rank correlation
with Benjamini-Hochberg FDR control, ranks hub nodes by degree and weighted
degree, and exports `network-edges.tsv`, `network-nodes.tsv`, and
`cooccurrence-network.png`.

## Key parameters

- `--min_prevalence`: minimum fraction of samples with non-zero abundance
  (default `0.2`) to avoid spurious zero-inflated correlations.
- `--top`: maximum number of highest-mean-abundance features retained when
  `--features` is not specified (default `100`).
- `--features`: optional explicit subset of features (e.g., the significant
  species returned by `mg_differential`), matching the biomarker co-occurrence
  workflow.
- `--min_abs_rho` / `--fdr`: edge significance cutoffs (defaults `0.4` and
  `0.05`).

## Interpretation and pitfalls

- Relative-abundance tables live on a simplex (sum to 1 or 100%), so Spearman
  correlations can induce spurious negative edges when one taxon expands
  dramatically. Treat negative correlations cautiously and validate critical
  hubs with compositionally aware tools (SparCC / SPIEC-EASI) if needed.
