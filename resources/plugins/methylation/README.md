# methylation · 0.1.0

Downstream DNA methylation analysis plugin (`toolPrefix: meth`, agent
`Methylation`, private locked Python 3.13 environment `plugin:methylation`) for
bisulfite/enzymatic sequencing (`nf-core/methylseq`, Bismark `.cov.gz`,
MethylDackel `.bedGraph`) and methylation arrays (Illumina 450K / EPIC
$\beta$-value matrices).

## Scope & division of work

- **Upstream alignment, deduplication, and methylation extraction**: run via
  `nf-core/methylseq`, `nf-core/bismark/*`, or `nf-core/methyldackel/*` wrappers
  (raw `.idat` fluorescence preprocessing for BeadChip arrays requires `minfi` /
  `sesame` externally before importing the normalized $\beta$ matrix).
- **Downstream statistical & biological chain (this plugin, 5 `meth_*` tools)**:
  1. `meth_inspect` — identify Bismark `.cov(.gz)`, MethylDackel `.bedGraph(.gz)`,
     or CpG/probe $\beta$ / $M$-value matrices, check coverage depth, bimodal
     $\beta$ range $[0, 1]$, and metadata alignment.
  2. `meth_merge` — merge per-sample Bismark `.cov(.gz)` or MethylDackel
     `.bedGraph(.gz)` files into coverage-filtered `beta_matrix.tsv` and
     `coverage_matrix.tsv` (`--min_cov`, `--min_samples_frac`).
  3. `meth_qc` — global methylation summary (hypo/intermediate/hypermethylated
     CpG fractions), $\beta$ density curves, logit $M$-value PCA, and global
     methylation shift / outlier diagnostics.
  4. `meth_dmr` — empirical-Bayes moderated linear modeling on logit $M$-values
     (with optional confounder covariates such as age, sex, batch, purity) paired
     with exact biological $\Delta\beta = \bar{\beta}_{\text{case}} - \bar{\beta}_{\text{control}}$
     filtering (`--delta_beta`, default `0.10`), plus coordinate-aware sliding
     merging of adjacent directional CpGs into **Differentially Methylated
     Regions (DMRs)** with Stouffer's Z regional significance and BH FDR.
  5. `meth_annotate` — annotate DMCs or DMRs against a reference GTF/GFF3
     (`Promoter`, `GeneBody`, `Downstream`, `Distal_Intergenic` + nearest TSS
     distance) and optional CpG Island BED (`Island`, `Shore`, `Shelf`,
     `OpenSea`), with optional GMT pathway ORA against the tested-CpG gene
     universe.
