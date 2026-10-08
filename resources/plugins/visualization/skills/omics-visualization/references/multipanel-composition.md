# Multi-panel composition

Use this reference only for labelled multi-panel figures or when several plots
must be rearranged into one scientific figure.

## Begin with the claim

Before choosing a grid, write a compact panel plan:

```text
Figure-level claim:
What could overturn the claim:
Panel sequence:
  a — scientific question | evidence role | necessary comparison
  b — scientific question | evidence role | necessary comparison
  c — scientific question | evidence role | necessary comparison
Hero panel:
Shared encodings:
```

Use one major claim per figure as a strong default. Panels may contain
subordinate findings only when they establish, compare, qualify, explain, or
bound that claim.

For Nature-style figures, the claim is normally carried by the caption and the
panel sequence, not by large in-panel titles. Keep panel-local text for axes,
direct labels, thresholds, units, and essential annotations.

## Give panels different jobs

Choose the smallest sufficient set. Useful evidence roles include:

| Role | Question |
|---|---|
| Overview or setup | What system, contrast, or population is being examined? |
| Primary evidence | Does the central effect or pattern exist? |
| Control or baseline | Does it exceed the relevant alternative? |
| Decomposition | Which components or strata account for the result? |
| Orthogonal validation | Does another assay or measure support it? |
| Boundary or failure case | Where does the conclusion weaken or stop? |
| Mechanistic evidence | What bounded explanation is supported? |

If removing a panel does not remove a unique inferential step, merge it, move
it to supplementary material, or omit it. Prefer evidence-role diversity over
showing the same result with several metrics.

## Assign visual hierarchy

- Give decisive evidence the hero position or largest area.
- Keep controls and validation panels quieter but readable.
- Preserve group colors, cohort names, denominators, uncertainty definitions,
  and comparator labels across panels.
- Share legends, axes, and color scales literally when panels use the same
  encoding. Do not duplicate small legends in every panel.
- Use equal panel sizes only when panels are peers.
- Use a schematic only when it is necessary to interpret the evidence.
- Panel letters express reading order, not universal panel functions. Use
  lowercase `a`, `b`, `c` by default for manuscript-like figures.

## Choose a layout

Read [layouts.yaml](layouts.yaml) and inspect only viable previews.

- `equal-2x2`: four peer panels.
- `hero-left-stack-right`: decisive evidence plus two supporting panels.
- `wide-top-two-bottom`: overview followed by two complementary details.

The layout source is an editable R composition template. Replace its demo panel
objects with the project-local plot objects and edit the visible Patchwork
design directly. Do not introduce a layout DSL.

Set the final physical size before adjusting text, legends, or gutters. Avoid
layout work that only looks balanced at an intermediate preview size.

## Final audit

Export at final physical size, run the layout's geometry audit, then follow
[rendered-layout-qa.md](rendered-layout-qa.md). An automated pass proves only
the declared geometry checks; it does not replace visual review.

This workflow is inspired by the claim-driven and rendered-geometry separation
in `nature-figure`, but it intentionally omits journal-specific panel counts,
legend limits, and submission rules.
