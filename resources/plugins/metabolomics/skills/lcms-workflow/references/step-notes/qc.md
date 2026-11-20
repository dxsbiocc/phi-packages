# mb_qc — QC-based feature filtering

## Principle

Pooled QC samples, injected repeatedly through the batch, measure technical
variation per feature. Features whose QC RSD is large are dominated by
instrument noise, not biology; features rarely detected in biological
samples carry little information. Filtering both is standard untargeted
practice (Dunn 2011; Broadhurst 2018 guidelines).

## Definitions used

- RSD = sample SD (ddof=1) / mean × 100 across QC samples. A feature needs
  ≥2 detected QC values; otherwise it fails. Mean ≤ 0 → fails.
- Detection = fraction of biological samples with value > 0 and non-NaN.
  QC and blank columns are excluded from this denominator.
- Defaults: RSD ≤ 30%, detection ≥ 0.8. Both are conventions, not laws.

## Interpretation

- Read `removed_by_rsd` vs `removed_by_detection`: RSD-dominated loss
  points at instrument/batch problems; detection-dominated loss at low
  abundance or poor extraction.
- The RSD histogram should pile up well below the threshold. A flat or
  bimodal distribution suggests mixed batches or mislabeled QCs.
- Typical untargeted plasma/urine datasets retain 50–80% of features.
  Losing most features is a red flag, not a result.

## Pitfalls

- No QC samples → the tool fails by design. Do not substitute blanks for
  QCs; a detection-only filter needs explicit user approval.
- Run before imputation/normalization: imputed values would fake stability.
- Blanks are reported, not subtracted; blank subtraction is separate manual
  work with its own assumptions.