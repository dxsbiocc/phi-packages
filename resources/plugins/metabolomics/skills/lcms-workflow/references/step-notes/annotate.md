# mb_annotate — m/z relationships and pseudospectra

## Principle

One compound ionizes into several features: adducts, isotopes, in-source
fragments. Deterministic mass-difference rules within ppm + RT tolerance
recover these relationships without any database (CAMERA's core idea,
Kuhl 2012 — simplified to m/z + RT, no peak-shape correlation).

## Rules implemented

- 13C isotope: Δ = 1.003355 Da.
- Adduct pairs (same neutral mass): pos [M+H]+ 1.007276, [M+Na]+ 22.989218,
  [M+K]+ 38.963158, [M+NH4]+ 18.033823; neg [M-H]- −1.007276,
  [M+Cl]- 34.969402, [M+FA-H]- 44.998201. A pair matches when
  Δmz ≈ delta(adduct2) − delta(adduct1).
- Neutral losses: H2O 18.010565, CO2 43.989829 (feature_2 is the higher
  mass, i.e. the parent).
- Tolerance: ppm relative to the observed higher m/z; RT window in table
  units (auto-converted from 0.2 min when the table is in seconds and no
  explicit `rt_tol` was passed).
- Pseudospectra = connected components over all matched relationships;
  consensus neutral mass comes from adduct-pair-derived neutral masses
  (median), falling back to base-adduct subtraction for isotope/loss-only
  groups.

## Interpretation

- This is **MSI level 3**: an annotated feature class, not an identification.
- Use groups to collapse redundancy: one compound ≈ one pseudospectrum,
  so pathway input should not count five adducts as five compounds.
- `ppm_error` near the tolerance edge is a weak match; tighten ppm and see
  if the relationship survives.

## Pitfalls

- Isobaric coincidences happen at wide ppm; RT proximity is the guard.
  Without an rt column the tool proceeds on m/z alone — weaker evidence,
  say so.
- Multiply charged ions, dimers ([2M+H]+), and custom adducts are not in
  the rule set; report them as gaps rather than forcing a match.
- Full CAMERA (correlation + peak shape across samples), formula
  assignment (Rdisop/enviPat), SIRIUS/MetFrag are gaps.