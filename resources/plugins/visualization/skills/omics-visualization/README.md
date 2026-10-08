**English** · [中文](README.zh-CN.md)

# Omics Visualization

A reusable agent skill for turning **existing omics results** into editable scientific figures. It is a catalog of canonical plot templates, not a statistics engine and not a place to invent charts from memory.

The skill exists so an agent can **choose a figure the way a careful analyst would**: state why the visualization helps, inspect the data distribution, name the claim, match the data shape, use real previews when they affect the decision, then adapt a known script instead of improvising ggplot2 from scratch. Its default polish follows a Nature-style philosophy: claim-first, restrained, compact, and readable at final size.

Agent operating rules live in [SKILL.md](SKILL.md). This page is the human-facing introduction.

## What it does

- Recommend a plot from a two-level catalog (family → template), using `use_when` / `avoid_when` rather than a filename guess.
- When the user only wants examples, show the installed `preview.png` files through `viz_examples`; no data simulation or rendering is needed.
- For an existing figure, edit its prepared script and rerender with the original inputs. For a user reference image, inspect it before making a new data figure.
- Use a lightweight contract router for preview-mode and common result-table shapes, so the agent can shortlist likely templates before opening full catalogs.
- Inspect the distribution and data contract before recommending a template.
- Show matching `preview.png` images when several candidates would change the scientific reading. Several candidates appear as a numbered list so you can compare them.
- Copy the chosen `plot.R` into **your project**, map columns, restyle, and render PDF / PNG / SVG.
- Apply named palettes from the OmicsAgent ColorPicker catalog instead of invented hex codes.
- Assemble labelled multi-panel figures when several plots must support one claim, then audit the layout.
- Move quickly when the catalog and data shape point to one clear template; show alternatives only when they change the scientific reading.

It does **not** run differential expression, enrichment tests, or causal
interpretation. If those results do not exist yet, the skill should say so
instead of fabricating numbers. A few templates compute **display** clustering
(`tree-dendrogram`, clustered heatmaps) from a supplied matrix; that is the
figure, not a substitute for an upstream analysis the user never ran.

## Design philosophy

### A figure is an argument

Titles, encodings, panel order, and uncertainty are part of the evidence, not decoration. The first question is *what claim should a reader be able to defend after seeing this figure?* The chart type follows that claim. A volcano plot is not the default for “omics data.”

### A plot serves the data

Visualization starts by stating what should become easier to see, compare, verify, or question. The agent inspects row counts, value ranges, category cardinality, sparsity, missingness, outliers, and matrix/network/hierarchy/time/genomic structure before recommending a template. If the requested chart would hide the important distribution, the skill should recommend a better-fitting catalog template and explain the tradeoff.

### Nature-style restraint

The baseline style is compact and evidence-led: minimal visual vocabulary, quiet secondary marks, explicit units and transformations, no decorative backgrounds, no rainbow palettes, and no labels that only work when zoomed. Journal-specific submission rules still come from the target venue; this skill supplies the design posture.

### Contract fast path, then catalog

Agents must not recommend a template from training data. For preview-mode and common shapes, they can first run the dependency-free router in [scripts/route_template.py](scripts/route_template.py), which uses generic entity/feature/category roles from [references/template_contracts.json](references/template_contracts.json). High-confidence results let the agent skip unrelated family catalogs, but the selected catalog entry and `plot.R` still get checked before rendering.

Contract edits are checked with [scripts/validate_template_contracts.py](scripts/validate_template_contracts.py), which verifies template ids, source paths, and preview paths against the catalogs. Its text report also shows contract coverage by family, so the remaining long tail can be expanded by priority instead of manually reading the whole catalog every time.

For low-confidence, novel, or publication-critical figures, agents read [references/plots.yaml](references/plots.yaml), open only that family’s catalog, and compare `preview`, `use_when`, `avoid_when`, and `input_shape`. Aliases (including Chinese names) help interpret a request; they are not enough to choose.

### Preview when it changes the decision

Use shipped previews to prevent imagined charts. If one template clearly matches the claim and data shape, record the id and proceed. If several viable templates would lead to different scientific readings, list up to four, each with `id`, title, one-line `use_when`, and the PNG.

Preview PNGs are **transparent**. A dark editor composites empty alpha as black. That black is the viewer, not the figure, and must not be “fixed.”

### Copy, then adapt

Bundled scripts are the canonical library. Ordinary use copies the selected `plot.R` into the project and edits that copy, while sourcing installed read-only helpers such as `scripts/lib/common.R` by absolute path. Copy `common.R` only for a one-off helper edit. Do not copy `references/` or palette catalogs into the project just to render; choose palette ids from the installed references. The installed skill is not the working file. There is no hidden plotting DSL: change `CONFIG` for columns and labels, `DATA PREPARATION` for reshape or order, and `PLOT` only when the geometry must change.

### One canonical implementation

The selected template decides the language. Current templates are R (`ggplot2` + `ggprism`). Do not rewrite a validated script into another language for stylistic uniformity, and do not add a second official implementation of the same `id`.

### Integrity over polish

Never invent sample sizes, tests, p-values, adjusted p-values, effect sizes, or biological meaning. Do not silently filter, impute, clip, or subsample. Report every material transform with before-and-after counts. Visual or statistical separation is not biological importance. Color-vision-safer defaults (`Qualitative.Safe`) beat decorative Artwork palettes unless you ask for them.

### Multi-panel is a claim, not a dashboard

A labelled figure starts from one major claim and a necessity test: if removing a panel removes no unique inferential step, merge it, move it, or drop it. Layouts come from [references/layouts.yaml](references/layouts.yaml). Geometry is audited; a failed audit blocks any claim that alignment passed. A single inset composite (circular tree plus outer enrichment track) is **one** template (`tree-enrichment-ring`), not a multi-panel layout.

## How selection works

```text
example request without data → viz_examples → shipped preview.png images → user choice

data + claim
    → purpose + distribution snapshot
    → route_template.py    (fast shortlist for preview/common shapes)
    → plots.yaml          (pick one family)
    → catalog/<family>.yaml
    → preview.png list    (only for material alternatives)
    → palettes.yaml       (named hex only)
    → copy plot.R locally
    → source installed helpers
    → nature-style pass   (final-size polish when relevant)
    → render, run qa_single_plot.py, and inspect the artifact
```

For a multi-panel request, write the figure-level claim and each panel’s evidence role first, then route every necessary panel through the same catalog.

## Library

About **149** single-plot templates in **12** families, plus a small layout set:

| Family | Typical use |
|---|---|
| `scatter` | Coordinates, embeddings, UMAP-in-circos rings, Voronoi tiles, volcano / zoom-inset volcano / MA, dumbbells, beeswarm, SVG glyphs, signed and circular lollipops, ternary compositions, pairs grids |
| `heatmap` | Matrices, clustering, cutree-block scores, correlation, association-dot / bubble / diagonal-split grids, Mantel-style displays, DE heatmap with enrichment zooms, term-by-group enrichment heatmap, saturation-mutation energy with mean ΔΔG overlay, gene × sample oncoprint, group-split circular heatmap |
| `bar` | Magnitudes, stacks, ranked enrichment bars and dotplots, category SVG marks, UpSet combination matrices, 2–4 set Venn circles, patient swimmer lanes |
| `boxplot` | Distributions, violins, rainclouds, paired contrasts |
| `line` | Time, rank, stacked area, survival, density ridges |
| `pie` | Part-to-whole, rings, Nightingale, circos doughnut with outer bars |
| `radar` | A few named numeric axes |
| `graph` | Node–link networks (supplied coordinates or `graph-stress`), weighted Circos chords |
| `sankey` | Conserved flows, including parallel-sets across categorical axes |
| `sunburst` | Hierarchy as rings whose **angle** encodes size |
| `tree` | Node-link trees, hclust dendrograms, treemap / pack / icicle, enrichment ring |
| `ideogram` | Gene models, chromosome G-banding, loci on chromosomes, density fill, genomic-locus coverage, protein mutation lollipops, dense-site dandelions, multi-track genomic Circos, two-genome synteny, genomic Circos heatmap, nested Circos zoom |

`circular` and `outer = "bar"` / `"point"` are CONFIG on the existing id, not extra templates. The SVG file on `scatter-svg` / `bar-svg-icon` is CONFIG too. `graph-force` keeps supplied x/y (`layout = "manual"`). Algorithm layout from an edge table is `graph-stress`.

Each template ships with `plot.R`, example data, and `preview.png`. Catalog entries record `id`, aliases, language, status, and the input shape the script expects.

Palettes live in the generated [palette index](references/palettes.yaml), with the full OmicsAgent set in the [color catalog](references/palettes/colors.json). These files are synced byte-for-byte from Phi's shared `resources/palettes/` catalog so the installed plugin stays self-contained; edit the shared source, not this copy. Defaults: `Qualitative.Safe` for groups, `Quantitative.BluGrn` for heatmaps, `Diverging.RdBu` for signed values, `Brand.Algolia` when matching the OmicsAgent product.

## Layout of the skill

```text
SKILL.md                         Agent contract
README.md / README.zh-CN.md      This introduction
references/
  template_contracts.json         Fast routing contracts for common shapes
  plots.yaml                     Family index
  palettes.yaml                  Generated recommended palettes
  palettes/colors.json           Generated full hex catalog
  catalog/<family>.yaml          Template records
  layouts.yaml                   Multi-panel layouts
  multipanel-composition.md
  nature-figure-principles.md
  rendered-layout-qa.md
scripts/
  route_template.py              Dependency-free preview router
  qa_single_plot.py              Lightweight single-plot artifact QA
  validate_template_contracts.py Contract/catalog consistency check
  lib/common.R                   Shared I/O and save helpers
  <family>/<template>/plot.R
  layouts/<layout>/compose.R
```

## Using it

In a compatible agent environment, invoke the **omics-visualization** skill when you have result tables or matrices and want a scientific figure. Bring the data (or a clear column contract). The agent should first state the display purpose, inspect the distribution, proceed when one template is clearly supported, and show previews for user selection only when viable candidates would tell materially different scientific stories.

To run a template yourself:

```bash
skill_run({ skill: "omics-visualization", script: "bar/basic/plot.R", args: ["scripts/bar/basic/example.tsv", "output.pdf"] })
```

Required R packages are listed in each script header. Missing packages should be reported, not installed silently.

## What this skill is not

- Not a replacement for DESeq2, limma, clusterProfiler, or any upstream test.
- Not a generic business-dashboard or AI-illustration kit.
- Not the workflow for authoring or mass-editing the template library itself.

[Switch to 中文](README.zh-CN.md)
