# Nature-Style Figure Principles

Use this reference when refining an omics figure for manuscript-like output or
when the user asks for Nature, Nature-like, or publication-ready polish. This
is a design philosophy, not a substitute for the target journal's current
submission instructions.

## Core Contract

- State the display purpose before plotting: what aspect of the data should be
  easier to see, compare, verify, or challenge?
- One figure, one defensible claim. A panel may qualify or bound the claim, but
  it should not start a second story.
- The visual form follows the evidence role. Choose geometry from the catalog
  before styling, and change styling before changing geometry.
- Clarity beats novelty. A familiar chart with precise encodings is better
  than a decorative chart that needs explanation.
- Final-size readability is the standard. A figure that only works zoomed in is
  not finished.

## Template Recommendation

Recommend from the data and purpose, not from the user's first chart noun.
Before committing to a template, check:

- whether the data are points, grouped distributions, ranked terms, a matrix,
  a hierarchy, a genomic interval set, a network, or a longitudinal record;
- whether numeric values are signed, bounded, sparse, zero-heavy, skewed, or
  dominated by outliers;
- whether categories are few enough for color, ordered enough for position, or
  too numerous for labels at final size;
- whether uncertainty, sample size, paired structure, censoring, or missingness
  must remain visible for the claim to be honest.

If several templates remain, show candidates only when they imply different
scientific readings or different overview/detail tradeoffs. Candidate notes
should name the template `id`, show the preview when available, and state both
the fit and the possible failure mode.

## Single-Plot Pass

Before final export, check the plot at its intended physical size:

- The main comparison is visible within three seconds.
- Axis labels include units or transformations, such as `log2 fold change` or
  `-log10 adjusted p`.
- Thresholds, confidence intervals, and uncertainty labels state exactly what
  they represent.
- Legends are short, non-redundant, and close enough to the data to avoid
  lookup work. Direct labels are preferred when they reduce clutter.
- Long gene, pathway, cohort, or cell-type names are wrapped, abbreviated with
  an explicit mapping, or moved to a more suitable orientation.
- Non-focal data are visible but quiet; focal data receive the strongest
  contrast.

## Multi-Panel Pass

For labelled figures, use [multipanel-composition.md](multipanel-composition.md)
first, then apply these finishing rules:

- Lowercase panel letters (`a`, `b`, `c`) are the default. Place them
  consistently and keep them separate from data labels.
- Shared legends, shared axes, and shared color scales should be shared
  literally, not recreated in every panel.
- Use equal panel sizes only for peer evidence. A hero panel should earn extra
  space by carrying the decisive comparison.
- Align plot areas, not just outer boxes. Axis text, color bars, dendrograms,
  and long labels should not create accidental panel hierarchy.
- Remove any panel whose result is already implied by another panel unless it
  provides a control, decomposition, validation, boundary, or mechanism.

## Efficient Application

- Apply this pass after template selection. Do not restart catalog routing
  unless the chosen template cannot support the scientific claim honestly.
- Prefer editing labels, scales, theme settings, legend guides, and output size
  over reworking data preparation.
- Rerender once per coherent set of visual edits, then inspect the artifact.
  Avoid tweak-render loops around changes that are invisible at final size.
- Keep the source editable and vector-first. Do not rasterize text, axes, or
  vector panels to solve layout problems.
