# Changelog — genomics

## 0.1.0 (2026-11-24)

Initial release: pure conda-forge downstream somatic cancer genomics plugin
(`toolPrefix: gen`).

- `gen_inspect`: GDC MAF validation, non-synonymous/silent breakdown, VAF/depth
  QC, clinical alignment, and optional CNV `.seg` inspection.
- `gen_summary`: cohort variant summary, 6-class substitution spectrum & Ti/Tv,
  MAF summary dashboard, and waterfall oncoplot (`Multi_Hit`).
- `gen_tmb_ith`: FoCR-compliant TMB (explicit `capture_size_mb` + `min_vaf`)
  and Mroz/Rocco `MATH` intratumor heterogeneity score with optional CNV segment
  exclusion and clinical group testing.
- `gen_interactions`: pairwise Fisher co-occurrence/mutual exclusivity (+ BH
  FDR), Sanchez-Vega (Cell 2018) 10 canonical oncogenic pathways, and
  cohort/clinical-group differential mutation testing (Fisher + OR + forest
  plot).
- `gen_signatures`: SBS96 trinucleotide context spectrum and `scipy.optimize.nnls`
  reference signature refitting with cosine similarity diagnostics.
- `gen_cnv`: CBS/CNVkit `.seg` Fraction Genome Altered (FGA), recurrent gain/loss
  chromosome summaries, and clinical group comparison.
