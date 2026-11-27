# Step note: TF regulon activity (`regulon`)

## Purpose
Infer single-cell transcription factor (TF) regulon activity per cell directly inside `plugin:scanpy` using either **AUCell** (recovery AUC among top expressed genes) or **ULM** (weighted univariate linear model $t$-statistic), stored in `adata.obsm["X_regulon"]`.

## How it works
- Reads unscaled log-normalized expression from `adata.layers["lognorm"]` (fallback `adata.X`).
- Matches TF target sets from a user `--regulons_path` (`source`/`tf`, `target`/`gene`, optional `weight`, or `.gmt`) or the curated 18-TF starter panel (`curated_tf_starter_v1`).
- Computes group-level differential regulon activity (two-sided Mann-Whitney U vs rest + BH FDR) across clusters/cell types when `group_by` exists in `adata.obs`.

## Reading the results
- `regulon_activity.tsv`: per-cell × TF activity matrix.
- `regulon_group_summary.tsv`: group means, cross-group z-scores, and BH FDR per `(group, tf)`.
- `figures/regulon-heatmap.png`: clustered group activity heatmap of top variable TFs.

## Pitfalls
- Expression-based regulons infer **transcriptional co-expression activity**, not physical DNA binding; confirm motif occupancy in chromatin accessibility (`single-cell-atac` or `chromatin-profiling`).
- With small starter target sets (< 5 expressed genes per TF), relax `regulon_min_targets` cautiously or pass a curated CollecTRI/DoRothEA table via `--regulons_path`.
