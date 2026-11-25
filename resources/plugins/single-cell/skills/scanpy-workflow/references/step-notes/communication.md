### Cell-cell communication (`communication`)

- **Purpose**: Infer putative paracrine/juxtacrine signaling axes (`sender -> receiver` via `ligand -> receptor`) across annotated single-cell populations.
- **Principle**: Uses the CellPhoneDB heteromeric subunit-minimum rule on the unscaled `lognorm` expression snapshot: every subunit of a multi-subunit ligand or receptor complex (`IL6R+IL6ST`, `TGFBR1+TGFBR2`, `ITGAV+ITGB1`) must reach `ccc_min_pct` expressed cells in its population, complex activity equals the limiting (minimum) subunit mean, and empirical specificity $p$-values are computed against `ccc_perms` cell-label permutations followed by BH FDR.
- **Interpretation**:
  - `communication-heatmap.png` shows which cell types act as dominant signaling hubs (`sender` rows vs `receiver` columns).
  - `communication-dotplot.png` ranks individual `ligand -> receptor` pairs by interaction strength (color) and permutation specificity (dot size $= -\log_{10}\text{FDR}$).
- **Pitfalls**:
  - Transcriptional co-expression measures **signaling potential**, not spatial proximity or protein phosphorylation; validate tissue contact in `spatial-transcriptomics` (`st_py_communication`).
  - When `interactions_path` is omitted, the step screens the built-in curated starter catalog (`curated_starter_v1`); supply a full OmniPath / CellPhoneDB TSV via `interactions_path` (or run `liana` in a project environment) for exhaustive receptorome discovery.
