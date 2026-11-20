# mb_preprocess — impute, normalize, transform

## Principle

Three deterministic operations in fixed order. Each is a modeling choice
the agent must make explicitly; defaults are all `none` so nothing happens
silently.

## Imputation (NaN and zero = missing)

- `half_min`: half the feature's minimum observed value — the common
  left-censored default (values below detection limit are small, not zero).
- `min`: the observed minimum; `zero`: literal zero (rarely appropriate
  before log).
- Features with no observed value stay missing and are counted in
  `n_missing_remaining`.
- Model-based left-censored methods (GSimp, QRILC) are a gap; heavy
  structured missingness deserves them.

## Normalization (per sample)

- `tic`: divide by column total. Simple; a few dominant features can drag
  every other feature down.
- `median`: divide by column median. Robust to dominant features.
- `pqn`: probabilistic quotient normalization (Dieterle 2006). Reference =
  sample with median TIC; per-sample factor = median of feature ratios to
  the reference (features where reference > 0). Standard for urine and
  other dilution-variable biofluids.
- `is`: divide by the internal-standard column (`is_col`), which is then
  dropped from the output. Requires a spiked IS measured per sample.
- QC-based drift correction (crmn, MetNormalizer, NormalizeMets) is a gap
  for long batches.

## Transform

- `log2` / `log10` stabilize variance and make fold changes additive.
  Refuses zeros/missing — impute first. After logging, `mb_univariate`
  needs `logged: true`.

## Interpretation

- `factors` is the audit trail: extreme factors (>5× or <0.2× the cohort
  median) mark suspect injections. Investigate before statistics.
- Keep the processing state in every report: "QC-filtered, half_min
  imputed, PQN normalized, log2" is one sentence that prevents most
  misreadings.