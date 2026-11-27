# Bulk RNA-seq Analysis Playbook

## 1. Upstream Wrapper Selection (`local/modules/differential-expression-*`)

| Scenario | Recommended wrapper | Notes |
| --- | --- | --- |
| Standard two-group or multi-factor raw counts (`n >= 3` per group) | `local/modules/differential-expression-deseq2` | Default choice; use `formula: "~ batch + condition"` when `rna_qc` shows a non-confounded batch effect |
| Large cohorts or strong dispersion / low-count emphasis | `local/modules/differential-expression-edger` | Quasi-likelihood F-test (`glmQLFTest`) on raw counts |
| Pre-normalized log-expression / microarray / heteroscedastic designs | `local/modules/differential-expression-limma` | Supports voom on counts or direct linear modeling on log-normalized matrices |
| Longitudinal / multi-timepoint trajectories | `local/modules/differential-expression-timeseries` | LRT / spline time-course modeling |

## 2. Downstream Call → Read → Decide Loop

1. **`rna_inspect`**: always call first. Confirm `inferred_kind == "raw_counts"` before DESeq2/edgeR, verify zero missing metadata samples, and inspect `batch_collinearity.confounded`.
2. **`rna_qc`**: inspect `pca_variance_explained` and `flagged_outliers`. If an outlier sample has low `mean_pearson_r` (< 0.80) due to failed library prep, confirm with the user before excluding it.
3. **Run DE wrapper** (`deseq2` / `edger` / `limma` / `timeseries`).
4. **Multi-contrast synthesis (`rna_compare`)**: when a study has ≥ 2 contrasts (e.g., DrugA vs Ctrl and DrugB vs Ctrl, or 6h vs 0h and 24h vs 0h), quantify `shared_up`, `shared_down`, and `discordant` genes rather than eyeballing separate volcano plots.
5. **Pathway interpretation (`rna_ora` + `rna_gsea` + `rna_pathway_score`)**:
   - Run `rna_ora` separately for `direction: up` and `direction: down` (mixing up and down in `both` dilutes directional pathways).
   - Run `rna_gsea` on the full unfiltered DE table ranked by `stat` to capture coordinated subtle shifts missed by hard DEG cutoffs.
   - Run `rna_pathway_score` when sample-level pathway heatmaps or covariate correlations are needed.
6. **TF/pathway regulon activity (`rna_regulon`)**:
   - Prefer `method: ulm` for independent pathway/TF testing and `method: mlm` when dissecting correlated regulons (e.g., CollecTRI/PROGENy) simultaneously.
   - Ensure gene symbols (not unmapped Ensembl IDs) are used; run `rna_id_map` first if needed.
7. **Differential alternative splicing (`rna_splicing`)**:
   - Pass all five rMATS `*.MATS.JC.txt` tables (`SE`, `A5SS`, `A3SS`, `MXE`, `RI`) or SUPPA2 `.dpsi` outputs together via `splicing_inputs` for unified FDR filtering, `frameshift_candidate` (`exon_len_bp % 3 != 0`) annotation, and multi-hit gene ranking.
