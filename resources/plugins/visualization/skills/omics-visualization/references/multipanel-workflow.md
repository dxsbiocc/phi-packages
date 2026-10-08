# Multi-panel composition workflow

Follow this for a labelled multi-panel figure, figure assembly, or rearrangement of
several plots. Code-span paths such as `references/plots.yaml` are relative to the
skill root; links on this page are relative to this file.

For a labelled multi-panel figure, figure assembly, or rearrangement of several
plots:

1. Read
   [references/multipanel-composition.md](multipanel-composition.md)
   and write the figure-level claim, panel sequence, evidence role of every
   panel, hero panel, and shared encodings before drawing.
2. Apply the necessity test. Merge, move, or omit a panel that adds no unique
   inferential step; do not arrange repeated metrics as a dashboard by default.
3. Route each necessary panel through `references/plots.yaml` and only the
   relevant family catalogs. Select one canonical template for each panel.
4. Read [references/layouts.yaml](layouts.yaml), inspect only viable
   layout previews, and choose a layout whose `use_when` matches the evidence
   hierarchy and whose `avoid_when` does not apply.
5. Read the selected `compose.R` and
   [scripts/layouts/lib/layout_common.R](../scripts/layouts/lib/layout_common.R)
   in full. Create a project-local copy under the active Phi project working
   directory containing the composition source and its selected panel sources;
   prefer sourcing installed read-only helpers by absolute path. Copy
   `common.R` or `layout_common.R` only when a one-off helper edit is required,
   and do not copy `references/` or palette catalogs.
6. Replace the layout template's demo panel objects with the real project-local
   plot objects. Edit the visible Patchwork design directly; do not create a
   layout DSL or rasterize vector panels merely to arrange them.
7. Render at final physical size. Run `audit_patchwork_layout()` and preserve
   its `.layout-audit.json` beside the figure. A `FAIL` or `NOT_AUDITABLE`
   result blocks any claim that alignment passed.
8. Follow
   [references/rendered-layout-qa.md](rendered-layout-qa.md), inspect
   the complete figure and every panel, then rerender and reaudit after any
   change to text, legends, axes, annotations, panel size, or layout.

Multi-panel delivery is complete only when every panel has a distinct evidence
role, the rendered geometry audit passes, and the final-size visual inspection
finds no unresolved clipping, collision, hierarchy, or reading-order defect.

