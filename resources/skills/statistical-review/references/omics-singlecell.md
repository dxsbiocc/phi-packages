# Omics, Single-Cell & Compositional Data Audit Branch

Load this reference ONLY when auditing bulk RNA-seq, single-cell/spatial transcriptomics, microbiome, or proteomics/metabolomics analyses.

## 1. Domain Redlines Matrix

| Flaw | Why It Breaks Inference | Required Fix |
| --- | --- | --- |
| **Cell-Level Pseudoreplication** | Running Wilcoxon/t-test across 40,000 pooled cells treats correlated cells from the same donor as independent (`P < 1e-300`) | Donor-level **pseudobulk** (`DESeq2`/`edgeR`/`limma-voom`) or GLMM with donor random effect |
| **Normalized Values in Count Models** | Feeding `TPM`/`FPKM`/`RPKM` or log-normalized matrices into `DESeq2`/`edgeR` breaks the negative-binomial variance model | Use **raw integer counts** (or `tximport` length-scaled counts) with size factors |
| **scRNA-seq Double Dipping** | Reporting `P`-values for cluster markers discovered by clustering on the exact same genes (**post-selection inference**) | Treat marker `P`-values as descriptive only, or validate on held-out split / independent cohort / `clusterpval` |
| **Compositional Closure Bias** | Standard t-test/Wilcoxon on cell-type proportions or microbiome relative abundances ignores sum-to-1 constraint | **CLR transform**, **ANCOM-BC**, **ALDEx2**, **scCODA**, or `propeller` / Dirichlet-multinomial |
| **Batch-Corrected Expression in DE** | Running DE on Harmony/ComBat/Seurat-integrated values distorts variance and inflates false positives | Run DE on **raw counts** with `batch` included as a covariate in the design matrix |
| **MNAR / LOD Left-Censoring** | Dropping or zero-imputing below-limit-of-detection proteomics/metabolomics values biases means and variance | Audit missingness mechanism; use left-censored imputation (`QRILC` / Tobit) + sensitivity analysis |

## 2. Worked Figure-Legend Upgrade

- *Before*: *"Fig. 2B. Expression of Gene X is higher in tumors (`*P < 0.05`, SEM shown)."*
- *After*: *"Fig. 2B. Normalized *GENE X* expression in matched tumor vs adjacent normal tissue (`n = 16` independent patients; dots = individual patients, center line = median, box = IQR, whiskers = `1.5 x IQR`). Two-sided paired Wilcoxon test (`W = 121`, `P = 0.0042`, `BH-FDR = 0.018`, median paired difference `+1.42 log2CPM`, `95% CI [0.58, 2.24]`)."*
