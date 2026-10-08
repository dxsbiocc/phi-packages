# Rendered layout QA

Use this after every multi-panel layout change and before delivery.

## Automated geometry gate

Each canonical layout calls `audit_patchwork_layout()` at the same width and
height used for export. It records final panel rectangles in points and checks
the relationships declared for that layout.

The audit covers:

- equal widths or heights for peer panels;
- shared left, right, top, and bottom boundaries;
- equal repeated horizontal or vertical gutters;
- outer alignment between a hero/overview panel and its supporting panels.

The default tolerance is 1.5 pt. Do not increase it to hide a single exception.
If a custom layout intentionally violates a canonical relationship, change the
layout or write a specific audit for that layout rather than claiming that a
canonical audit passed.

Audit statuses:

- `PASS`: every declared relationship is within tolerance.
- `FAIL`: at least one reliable geometric relationship is outside tolerance;
  fix and rerender before delivery.
- `NOT_AUDITABLE`: panel viewports could not be measured; do not claim that
  alignment passed.

## Final-size quality gate

Judge the exported figure at the physical size where it will be read. For
Nature-style or manuscript-like output, confirm:

- tick labels are readable without zooming, with axis and legend text at least
  as prominent as the smallest important data label;
- thin lines, grid lines, borders, and dendrogram branches survive reduction
  without overpowering the data;
- panel letters are lowercase, bold, consistently placed, and separate from
  data labels;
- units, transformations, p-value type, uncertainty definition, and group
  denominators are explicit where they affect interpretation;
- color encodings remain interpretable in grayscale or with redundant shape,
  linetype, label, ordering, or annotation.

## Visual inspection

At final physical size, inspect the complete figure and every panel:

- reading order and panel letters;
- whether the strongest evidence is visually dominant;
- axis and plot-area alignment;
- repeated gaps and outer margins;
- title, label, annotation, and legend collisions;
- clipped text, data, or significance marks;
- inconsistent fonts, colors, group mappings, or uncertainty encodings;
- redundant legends and axes;
- decorative color, background, or annotation that does not carry evidence;
- heatmap legends, network labels, and long category names that distort the
  intended grid;
- whether each panel remains interpretable after the final reduction.

Any visual change to text, legends, axes, annotations, panel size, or layout
invalidates the prior review. Rerender, rerun the geometry audit, and inspect
again.
