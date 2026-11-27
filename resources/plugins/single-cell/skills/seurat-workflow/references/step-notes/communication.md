# Step note: Cell-cell communication (`communication`)

## Purpose
Screen heteromeric ligand-receptor signaling potential across Seurat clusters or cell types using pure R (`SeuratObject::LayerData`, `Matrix`, `ggplot2`), mirroring the Python `sc_py_communication` CellPhoneDB subunit-minimum permutation model without external LIANA dependencies.

## How it works
- Computes sender ligand and receiver receptor subunit-minimum expression and detection fractions from the unscaled normalized `data` layer.
- Gates pairs on `ccc_min_pct` (default `0.10`), permutes cell-type labels `ccc_perms` times (default `200`), and adjusts empirical p-values by Benjamini-Hochberg FDR (`ccc_fdr`, default `0.05`).

## Reading the results
- `interactions.tsv`: per `(sender, receiver, ligand, receptor)` mean score, detection percentages, empirical p-value, and BH FDR.
- `communication-heatmap.png`: significant interaction count matrix across sender and receiver populations.
- `communication-dotplot.png`: top ranked signaling pairs colored by `-log10(FDR)`.

## Pitfalls
- mRNA co-expression reflects signaling potential, not physical contact (validate spatially via `spatial-transcriptomics`).
- For multi-method consensus (`CellChat` / `NicheNet` / `LIANA` OmniPath), run `sc_env_probe` and request a project environment.
