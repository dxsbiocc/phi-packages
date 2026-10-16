---
name: create-template
description: >-
  Create, extend, register, and validate reusable omics figure templates in the
  sibling omics-visualization skill (plot.R, example tables, catalog, and
  routing contracts). Use when the user attaches a figure and asks to reproduce
  it as a template, or says 做成模板, 复现为绘图模板, 保存模板, 扩展模板,
  新增模板, or 复现模板, wants a new template-id for a geometry CONFIG flags
  cannot express, or asks to promote a one-off plot into the catalog (template
  authoring, register template). Do not use for ordinary plotting of user data
  with an existing template, revisions, example previews, or reference-guided
  user-data plots; those are omics-visualization.
phi:
  attachTo: [Visualization]
---

# Create Template

Create and register reusable plot templates in the sibling
[`omics-visualization`](../omics-visualization/SKILL.md) skill (same
visualization plugin: `skills/omics-visualization/`). That directory is the
installed omics-visualization skill root. Rendering user data with an existing
template is not this skill; that is `skill://omics-visualization`.

## When to use

Trigger on any of:

1. **Image → template**: the user attaches a figure and asks to reproduce it as
   a template / 做成模板 / 复现为绘图模板.
2. **Save / register**: the user says 保存模板, 新增模板, 扩展模板, 复现模板,
   or asks to promote a working plot into the catalog.
3. **Geometry gap**: the closest catalog entry cannot express the claim without
   a new visual grammar (not just a layout flag).

Do **not** trigger for: pick an existing template and plot a user table; polish
labels on a project-local copy; show example previews; use a reference image to
style a plot of the user's data; upstream statistics. Those stay on
`omics-visualization`.

## Where files go

Edit the sibling skill `../omics-visualization/` only. Catalog YAML, routing
contracts, `plot.R`, example tables, and previews live there. This skill has no
`scripts/` directory. Do not copy `validate_template_contracts.py`,
`qa_single_plot.py`, or template `plot.R` files into `create-template`.

Resolve the installed create-template directory (the skill you are reading) and
edit its sibling. Do not guess a source-repository checkout, and do not treat a
project-local `plot.R` as registered.

## Decision: extend vs new id

Read `../omics-visualization/references/plots.yaml` and the relevant
`../omics-visualization/references/catalog/<family>.yaml` before inventing an id.

| Choice | When |
|--------|------|
| **CONFIG / layout flag** | Same claim and geometry; only glyph, circular/polar, palette, or column names change. Not a new id. |
| **Extend sibling script** | Same template-id; schema honesty fix (e.g. require `pair_id` instead of row-order pairing). |
| **New template-id** | Materially different claim or geometry (e.g. dual-ring hub network ≠ single-ring `graph-circular`). |
| **Refuse / redirect** | Requested chart is decorative, not data-backed, or belongs to another family already covered. |

Never add a sibling id that only changes `layout`, `circular`, force-algorithm
alias (`kk` / `fr` / `lgl`), or SVG glyph file.

## Workflow

Copy this checklist and track it:

```
Author progress:
- [ ] 1. Decode target (image or request) → claim + encodings + input shape
- [ ] 2. Choose family and extend-vs-new-id
- [ ] 3. Implement ../omics-visualization/scripts/<family>/<slug>/plot.R (+ example tables)
- [ ] 4. Register catalog, contracts, aliases, mis-routes in that sibling skill
- [ ] 5. Render preview; QA; validate contracts from the omics-visualization skill root
- [ ] 6. Fix honesty / palette / legend issues; re-render
```

### 1. Decode the target

If the user attached an image, **read the image file** (do not rely only on a
text description). Extract:

- One-sentence **figure claim**
- Marks: nodes/points/bars/ribbons/boxes/…
- Encodings: colour / size / linewidth / shape / position / facets
- Required **input tables** and columns (observation unit first)
- What must be **supplied** vs never computed in the script

Write a short data contract before coding. For image reproduction, invent
minimal example data that exercises every encoding; do not invent p-values or
edges the figure does not support.

Details: [references/image-decode.md](references/image-decode.md).

### 2. Mirror a sibling template

Pick the closest existing `../omics-visualization/scripts/<family>/*/plot.R`
and match its structure:

1. Banner: `Template-ID`, purpose, inputs, dependencies, example CLI,
   scientific assumptions
2. Locate `lib/common.R` via `--file=` walk-up
3. `parse_io_args(...)` then top-level `config <- list(...)`
4. `load_packages(...)` — fail if missing; never auto-install
5. DATA PREPARATION: `require_columns`, type checks, early `stop()` on bad
   contracts; display-only transforms only
6. PLOT: ggplot2 / ggraph / ComplexHeatmap / circlize as the family uses
7. SAVE: `save_ggplot(...)` (or family equivalent)

CLI shape examples, run from the template directory or with paths under the
omics-visualization skill root:

```bash
Rscript plot.R example.tsv output.pdf
Rscript plot.R nodes.tsv links.tsv output.pdf
```

Rules:

- User-editable mappings live in `CONFIG`; do not invent a JSON plot DSL.
- No upstream DE / enrichment / community detection / correlation inference
  unless the sibling already documents that display transform.
- Prefer explicit ids (`pair_id`, `ring`, `axis`) over row-order pairing.
- Colours: `palette_colors("Qualitative.Safe"|"Qualitative.Bold"|...)` or hex
  copied from `../omics-visualization/references/palettes.yaml` /
  `../omics-visualization/references/palettes/colors.json`.
  **Do not invent or interpolate hex.**

### 3. Example data + preview

Ship next to `plot.R` inside the sibling skill:

- `example.tsv` / `nodes.tsv`+`links.tsv` / family sidecars
- Headers matching default `config$columns`
- Enough rows to show every legend class and encoding

Render into that template directory:

```bash
Rscript plot.R <inputs...> preview.png
Rscript plot.R <inputs...> preview.pdf   # optional but preferred
```

### 4. Register

Minimum registration set. Paths are inside the sibling skill
`../omics-visualization/` (the installed omics-visualization skill root):

1. `references/catalog/<family>.yaml` — new entry: `id`, `title`, `aliases`,
   `source`, `examples`, `preview`, `use_when`, `avoid_when`, `input_shape`,
   `language`, `status`
2. Update nearby siblings' `avoid_when` so agents do not mis-route
3. `references/template_contracts.json` — same `id`, `preferred_shapes`,
   `roles` (only known role names from `scripts/route_template.py` in that
   skill), `intent_keywords` (include Chinese aliases), `use_when`,
   `avoid_when`
4. `references/plots.yaml` — family aliases if users will say a new Chinese name
5. `SKILL.md` **Easy mis-routes** — one line when the new id is easy to confuse

Id form: `<family>-<slug>` with underscores in the directory,
hyphens in the id (`graph/circular_concentric` → `graph-circular-concentric`).

Role names must exist in the contract validator (e.g. `pair_id`, `category`,
`numeric_y`, `weight`, `magnitude`, `axis_entity`, `source_entity`,
`target_entity`). Unknown roles fail validation.

Full checklist: [references/checklist.md](references/checklist.md).

### 5. Validate

Run these from the installed omics-visualization skill root (sibling
`../omics-visualization/`). They are that skill's scripts. Do not create them
under `create-template`.

```bash
python3 scripts/validate_template_contracts.py --json
python3 scripts/qa_single_plot.py scripts/<family>/<slug>/preview.png --json
```

Contracts must report `"ok": true` with no errors. Re-read the preview image
and fix: clipped legends, dishonest encodings, missing required columns,
palette policy violations, transparent-on-dark illegibility when the reference
figure is white-background.

### 6. Delivery

Tell the user:

- template id
- input schema (columns / sidecars)
- path to `plot.R` and preview inside the omics-visualization skill
- how to run: `Rscript plot.R ... output.pdf`
- extend vs new-id rationale (one sentence)

## Anti-patterns

- One-off project `plot.R` without catalog/contract registration when the user
  asked to save a template
- New id for a CONFIG-only variant
- Pairing or grouping inferred from sort order when an id column exists or
  should exist
- Invented hex, rainbow scales, or computing statistics the template claims
  are “supplied”
- Skipping `python3 scripts/validate_template_contracts.py` after editing
  contracts (run it from the omics-visualization skill root)
- Adding scripts under `create-template`
- Editing the installed omics-visualization skill during ordinary user-data
  plotting. Authoring is the exception, and it edits that sibling library.
  Ordinary plotting uses `viz_*` tools and a project-local copy.

## Relationship to omics-visualization

| Skill | Job |
|-------|-----|
| `omics-visualization` (`skill://omics-visualization`) | Choose an existing template, copy it into the project, and plot user data |
| `create-template` (`skill://create-template`) | Add or extend catalog templates, example data, contracts, and previews in the installed omics-visualization library |
