# Step note: GC-matched chromVAR TF motif deviations (`tf-deviation`)

## Purpose
Quantify per-cell transcription factor motif accessibility deviations ($z$-scores) using GC-content and mean-accessibility 5×5 quantile-matched background peak resampling (Schep et al. 2017 chromVAR model), cleanly implemented in Python without R/ArchR dependencies.

## How it works
- Computes expected per-cell, per-peak accessibility $E_{ij} = n_i e_j$.
- Bins peaks into $5 \times 5$ `(GC, log10(mean_access))` quantiles and draws `n_background` matched background sets per TF.
- Subtracting background deviation mean corrects technical GC/depth bias; dividing by background SD produces per-cell deviation $z$-scores (`obsm["X_tf_deviation"]`) and across-cell variability SD.

## Reading the results
- `tf_deviations.tsv`: cell × TF deviation $z$-score matrix.
- `tf_variability.tsv`: TFs ranked by cell-to-cell deviation SD and `mean_abs_z`.
- `figures/tf-deviation-summary.png`: top variable TF barplot + per-cluster mean z-score heatmap.

## Pitfalls
- Motif deviation reflects **accessibility variation at motif occurrences**, not direct physical TF ChIP occupancy.
- Supply `peak_gc` or `genome_fasta` whenever possible so GC matching is active rather than accessibility-only fallback.
