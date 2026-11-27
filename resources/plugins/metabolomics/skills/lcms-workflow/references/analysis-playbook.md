# LC-MS metabolomics playbook

The decision loop for every step: **call → read the JSON → decide → check red
lines**. Tools are deterministic; you (the agent) choose steps and thresholds.
This chain starts at an exported feature table, not at raw instrument files.

## 0. Route check

- Question is "which features differ / what are they / which pathways light
  up" on an untargeted LC-MS feature table → this skill.
- Raw files (`.raw`, `.mzML`, `.d`) → upstream peak picking is required and
  is **not** in this plugin (see §1). Do not pretend to read them.
- GC-MS feature tables work mechanically (mz/rt columns), but retention
  indices — the GC-MS identification anchor — are not supported. Say so.
- NMR spectra, targeted MRM quantification, and lipid-species naming
  (e.g. PC(16:0/18:1)) are out of scope. Report them as gaps.

## 1. Upstream (a gap, not a wrapper)

msconvert (ProteoWizard), xcms (centWave / matchedFilter), MZmine, and
MS-DIAL perform conversion, peak picking, alignment, and gap filling. None
has a wrapper in this repository today. Ask the user for:

1. the exported feature table (xcms `featureValues` + `featureDefinitions`,
   MZmine feature list, MS-DIAL alignment result, or MetaboAnalyst table),
2. sample metadata with group labels,
3. which samples are QC / pool / blank (names usually carry `QC`/`pool`),
4. optionally: MS2 spectra (MGF), a spectral library (MSP/MGF), a
   pathway→compound mapping, a compound mass database.

Confirm the RT unit at this point: xcms exports seconds, MZmine/MS-DIAL
minutes. `mb_inspect` infers it from the median RT; verify against the
vendor export before trusting the inference.

## 2. Inspect — `mb_inspect`

- Call on the file or directory. Read `files[].kind`, `format`, `rt_unit`,
  `qc_samples`, `missing_rate`, `ready_steps`.
- `metaboanalyst` format means samples-as-rows: transpose before anything
  else (the tool reports it; it does not transpose for you).
- `unknown` → open the file yourself (read tool) before guessing.
- Red line: never guess which samples are QC. If inspect finds none, §3
  cannot run as designed.

## 3. QC filtering — `mb_qc`

- Convention (Dunn 2011 / Broadhurst 2018 practice): keep features with
  QC RSD ≤ 30% and detection rate ≥ 0.8 across biological samples. These
  are starting points; tighten (20%) for large cohorts, loosen consciously
  and report the choice.
- Read `removed_by_rsd` / `removed_by_detection` and the RSD histogram.
  Losing >50% of features usually means a mislabeled QC set or a failed
  batch — stop and investigate, do not just relax thresholds.
- No QC samples → the tool fails. Report the gap; a detection-rate-only
  filter is a conscious fallback the user must approve, not a default.
- Blank samples are reported, not subtracted. Blank subtraction is a
  decision with its own assumptions; flag it as manual work.
- Run on **raw intensities**, before imputation or normalization.

## 4. Preprocessing — `mb_preprocess`

- Fixed order: impute → normalize → transform. NaN and zero both count as
  missing (below detection limit).
- Impute: `half_min` is the common left-censored default; `min` and `zero`
  are alternatives. Model-based left-censored imputation (GSimp, QRILC) is
  a gap — say so when missingness is heavy and structured.
- Normalize: `tic` (total intensity) is simple but sensitive to a few huge
  features; `median` is robust; `pqn` (median-of-ratios to a reference
  sample, Dieterle 2006) is the usual choice for urine/biofluids; `is`
  divides by an internal-standard column when one was spiked in.
  For long multi-batch runs showing injection-order signal drift, run
  `mb_drift_correct` (QC-RLSC robust LOESS) on raw intensities before
  `mb_qc` / `mb_preprocess`; internal-standard-only `crmn` normalization
  remains a gap.
- Transform: `log2` for statistics and figures. The tool refuses log on
  zeros/missing — impute first. After a log transform, pass `logged: true`
  to `mb_univariate`.
- Read `factors`: a sample whose normalization factor is an extreme
  outlier (>5× or <0.2× the cohort median) is a suspect injection — flag
  it to the user before statistics.

## 5. Multivariate — `mb_multivariate`

- Look at PCA first (unsupervised): it shows batch effects and outliers
  without using labels. PLS-DA is supervised and **will** separate any
  labels given enough components — separation alone proves nothing.
- Scaling: `uv` (unit variance) gives every feature equal weight;
  `pareto` keeps some intensity weighting (common in metabolomics);
  `none` lets high-abundance features dominate. State the choice.
- Trust PLS-DA only with the permutation result: `q2_p` and `r2y_p` from
  seeded label permutations. Q2 ≤ 0 or a non-significant permutation p
  means the model is overfit — report that honestly, never hide it.
- VIP > 1 is the conventional importance cut; combine with univariate FDR
  and fold change, never VIP alone.
- Small n (<10/group) makes Q2 unstable; the tool uses LOO CV (7-fold
  above n=50). Say which scheme ran.
- OPLS-DA (ropls / mixOmics) is a gap — do not fake it with PLS-DA output.

## 6. Univariate — `mb_univariate`

- Two groups: `wilcoxon` (robust, default) or `welch` (more power when
  roughly normal, e.g. after log). More than two: `kruskal`.
- Significance = BH FDR < 0.05 **and** |log2FC| ≥ 1 by convention; both
  thresholds are decisions, pass them explicitly.
- `logged: true` after a log transform (FC = 2^mean difference). Without
  it, FC is a raw mean ratio — do not mix the two interpretations.
- With n=3/group, Wilcoxon cannot reach p<0.05; a null result there is a
  power statement, not evidence of no effect. Say so.
- Run on biological samples only — drop QC/blank columns first.

## 7. Annotation — `mb_annotate`

- Deterministic m/z relationships: 13C isotopes, adduct pairs, neutral
  losses, grouped into pseudospectra by shared neutral mass + RT.
- Tolerances: 10 ppm and 0.2 min are starting points for Orbitrap/TOF;
  widen for low-resolution instruments, tighten for high. The RT window
  auto-converts when the table is in seconds — verify `rt_unit` in the JSON.
- This is **MSI level 3** (annotated feature class), not identification.
  CAMERA-style correlation + peak-shape grouping, Rdisop/enviPat formula
  assignment, and SIRIUS/MetFrag are gaps.
- Use the pseudospectra to collapse redundancy before pathway analysis:
  one compound should not count as five adduct features.

## 8. MS2 matching — `mb_ms2_match`

- Library provenance is part of the result: record which MSP/MGF (MassBank,
  MoNA, GNPS, NIST, in-house) and its version/license. Confirm with the
  user before using a commercial library.
- Cosine ≥ 0.7 with ≥ 3–5 matched fragments is a common reporting floor;
  the defaults are starting points. Precursor 10 ppm, fragments 0.02 Da
  suit high-resolution MS2.
- A library hit is **MSI level 2** (putative). Level 1 requires an
  authentic standard run on the same method (RT + MS2 match). Never
  upgrade a level-2 hit to "identified" in prose.
- No MS2 in the data → this step is impossible; say so instead of
  guessing from m/z alone.

## 9. Pathway enrichment — `mb_pathway`

- No built-in database. The user supplies the pathway→compound mapping
  (KEGG, Reactome, HMDB, MetaboAnalyst libraries) and, for `masses` mode,
  a compound mass table. Confirm licensing and version; cite them.
- `ids`: use when features are already putatively identified (level 2/3
  with names). Background = all compounds in the mapping — keep it
  consistent with what the platform could have detected.
- `masses`: mummichog-style (Li 2013) when most features are unidentified:
  m/z → adduct hypotheses → candidate compound sets, permutation over
  random m/z lists from the **full detected background**. The background
  must be all features, not just the significant ones.
- Enrichment is a hypothesis generator on ambiguous matches; report
  overlap counts and empirical p, and keep the MSI caveat attached.
- FELLA-style network diffusion and topology-aware enrichment are gaps.

## 10. Figures — `mb_plot`

- Heatmap: top N by variance, log10 + clustering for intensity data.
  Clustering on raw intensities mostly ranks abundance — pair with log10.
- Score plots come from `mb_multivariate`; volcano from `mb_univariate`.
  Every figure has a `.phi-artifact.json` descriptor — cite the artifact,
  do not regenerate numbers by eye.

## Red lines (report, never fake)

- Upstream peak picking / alignment (msconvert, xcms, MZmine, MS-DIAL).
- OPLS-DA; multi-batch crmn modeling; GSimp/QRILC imputation.
- CAMERA full pseudospectra, formula assignment (Rdisop/enviPat),
  SIRIUS / MetFrag / GNPS networking.
- GC-MS retention indices; NMR; targeted MRM; lipid-species naming.
- Level-1 identification without a standard; pathway results without a
  user-supplied, licensed mapping.

## Reporting rules

- Every JSON return is the audit record: cite its counts and thresholds.
- State the processing state of every matrix (raw / QC-filtered /
  imputed+normalized+logged) when quoting a number.
- Attach the MSI level to every identification claim.
- List gaps explicitly at the end of the analysis, not buried in prose.