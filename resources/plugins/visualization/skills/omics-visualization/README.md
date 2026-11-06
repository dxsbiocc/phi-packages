**English** · [中文](README.zh-CN.md)

# Omics Visualization

A reusable agent skill for turning **existing omics results** into editable scientific figures. It is a catalog of canonical plot templates, not a statistics engine and not a place to invent charts from memory.

The skill exists so an agent can **choose a figure the way a careful analyst would**: state why the visualization helps, inspect the data distribution, name the claim, match the data shape, use real previews when they affect the decision, then adapt a known script instead of improvising ggplot2 from scratch. Its default polish follows a Nature-style philosophy: claim-first, restrained, compact, and readable at final size.

Agent operating rules live in [SKILL.md](SKILL.md). This page is the human-facing introduction.

## What it does

- Recommend a plot from a two-level catalog (family → template), using `use_when` / `avoid_when` rather than a filename guess. The agent chooses the template; the tools look up, count, copy, and run.
- When the user only wants examples, look up installed `preview.png` files with `viz_catalog` by alias or family; no data simulation or rendering is needed.
- For an existing figure, edit its prepared script and rerender with the original inputs. For a user reference image, inspect it before making a new data figure.
- Count the table with `viz_profile` before choosing: column types, missing values, ranges, and identifier overlap.
- Look up one family or an exact name with `viz_catalog`. It never ranks or guesses.
- Show matching `preview.png` images when several candidates would change the scientific reading. Several candidates appear as a numbered list so you can compare them.
- Copy the chosen `plot.R` into **your project** with `viz_prepare`, map columns, restyle, and render PDF / PNG / SVG with `viz_render`.
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

### Agent chooses, tools look up

Agents must not recommend a template from training data, and scripts must not rank templates from column names or keywords. The agent reads the purpose, the `viz_profile` facts, and one family from `viz_catalog`, then compares `use_when`, `avoid_when`, and `input_shape`. `viz_catalog` looks up by family, exact id, or exact alias (including Chinese names). It does not score or shortlist.

[references/template_contracts.json](references/template_contracts.json) lists the input roles each template expects. `viz_prepare` checks an agent-supplied `role=column` mapping against those roles; it does not choose the template. Catalog edits are checked with [scripts/validate_template_contracts.py](scripts/validate_template_contracts.py), which verifies template ids, source paths, preview paths, and role names against the catalogs.

For a named chart, look it up with `alias`. For a new figure, read [references/plots.yaml](references/plots.yaml), open only that family’s catalog, and choose. An alias helps interpret a request; it is not enough to choose when several templates would change the reading.

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
example request without data
    → viz_catalog (alias or family)
    → shipped preview.png images
    → user choice

data + claim
    → purpose
    → viz_profile          (count the table)
    → plots.yaml           (pick one family)
    → viz_catalog family   (read use_when / avoid_when / input_shape)
    → preview.png list     (only for material alternatives)
    → palettes.yaml        (named hex only)
    → viz_prepare          (copy plot.R; check role=column)
    → source installed helpers
    → nature-style pass    (final-size polish when relevant)
    → viz_render, inspect the artifact
```

For a multi-panel request, write the figure-level claim and each panel’s evidence role first, then route every necessary panel through the same catalog.

## Library

**211** single-plot templates in **12** families, plus a small layout set. Some templates also ship `preview_<name>.png` CONFIG variants:

| Family | Typical use |
|---|---|
| `scatter` | Coordinates, embeddings, UMAP atlas / circos / marker-gene dots, Voronoi tiles, volcano / multi / barcode / GO-shoulder / MA, Manhattan / QQ, 1D hist / density, 2D density, residuals, xy-error, zoom insets, dumbbells, beeswarm (no boxes), SVG glyphs, lollipops, ternary, pairs, grouped fits, pathway facets, Zi–Pi, variate panels |
| `heatmap` | Matrices, clustering, cutree-block scores, correlation, association-dot / bubble / diagonal-split / observation-dot grids, Mantel-style displays, DE heatmap with enrichment zooms, term-by-group enrichment, Circos annotation / cohort / split rings, alteration-frequency cohort, saturation-mutation energy, oncoprint |
| `bar` | Magnitudes, stacks, stack+doughnut, ranked enrichment bars / cluster / lollipop, cell-communication bars, radial donut / line, multi-axis line+bar, Circos links, category SVG marks, UpSet, Venn, flower / petal sets, swimmer lanes |
| `boxplot` | Distributions, outline-box beeswarm / notched / cohort boxes, violins, rainclouds, letter-value boxes, paired and butterfly contrasts |
| `line` | Time, rank, stacked area, survival, density ridges, GSEA, difference ribbons, ROC, spiral plots |
| `pie` | Part-to-whole, rings, Nightingale / petal pie, circos doughnut with outer bars |
| `radar` | A few named numeric axes, including significance-marked profiles |
| `graph` | Node–link networks (supplied coordinates or `graph-stress`), TF–target regulons, signed multi-omics networks, module co-occurrence panels, weighted / bipartite Circos chords |
| `sankey` | Conserved flows, parallel-sets, equal-width protein–trait alluvial maps, enrichment flow with a ranked bar |
| `sunburst` | Hierarchy as rings whose **angle** encodes size, including petal taxonomies |
| `tree` | Node-link trees, hclust dendrograms, treemap / pack / icicle, enrichment ring, circular cell-annotation tree, daylight and phylo-ring phylogenies |
| `ideogram` | Gene models, chromosome G-banding, loci / multi-trait loci, density fill, genomic-locus coverage, protein mutation lollipops, dandelions, Circos / synteny / linear riparian / subgenome Circos, nested Circos zoom |

`circular` and `outer = "bar"` / `"point"` are CONFIG on the existing id, not extra templates. The SVG file on `scatter-svg` / `bar-svg-icon` is CONFIG too. `graph-force` keeps supplied x/y (`layout = "manual"`). Algorithm layout from an edge table is `graph-stress`.

Each template ships with `plot.R`, example data, and `preview.png`. A CONFIG option that changes the figure also ships as `preview_<name>.png` and is listed under `variants` in the catalog. `viz_catalog` `ids` returns those images as extra `preview_markdown`. Catalog entries record `id`, aliases, language, status, and the input shape the script expects.

Palettes live in the generated [palette index](references/palettes.yaml), with the full OmicsAgent set in the [color catalog](references/palettes/colors.json). These files are synced byte-for-byte from Phi's shared `resources/palettes/` catalog so the installed plugin stays self-contained; edit the shared source, not this copy. Defaults: `Qualitative.Safe` for groups, `Quantitative.BluGrn` for heatmaps, `Diverging.RdBu` for signed values, `Brand.Algolia` when matching the OmicsAgent product.

## Layout of the skill

```text
SKILL.md                         Agent contract
README.md / README.zh-CN.md      This introduction
references/
  template_contracts.json        Input roles each template expects
  plots.yaml                     Family index
  palettes.yaml                  Generated recommended palettes
  palettes/colors.json           Generated full hex catalog
  catalog/<family>.yaml          Template records
  layouts.yaml                   Multi-panel layouts
  multipanel-composition.md
  nature-figure-principles.md
  rendered-layout-qa.md
scripts/
  viz.py                         catalog / profile / prepare / render
  viz_catalog.py                 Exact catalog lookup
  viz_profile.py                 Table counts and column checks
  viz_templates.py               Copy and parse plot.R
  qa_single_plot.py              Lightweight single-plot artifact QA
  validate_template_contracts.py Catalog and role consistency check
  lib/common.R                   Shared I/O and save helpers
  <family>/<template>/plot.R
  layouts/<layout>/compose.R
```

## Using it

In a compatible agent environment, invoke the **omics-visualization** skill when you have result tables or matrices and want a scientific figure. Bring the data (or a clear column contract). The agent should first state the display purpose, count the table with `viz_profile`, look up one family with `viz_catalog`, proceed when one template is clearly supported, and show previews for user selection only when viable candidates would tell materially different scientific stories.

To run a template yourself:

```bash
skill_run({ skill: "omics-visualization", script: "bar/basic/plot.R", args: ["scripts/bar/basic/example.tsv", "output.pdf"] })
```

Required R packages are listed in each script header. Missing packages should be reported, not installed silently.

## What this skill is not

- Not a replacement for DESeq2, limma, clusterProfiler, or any upstream test.
- Not a generic business-dashboard or AI-illustration kit.
- Not the workflow for authoring or mass-editing the template library itself. That is the sibling skill `create-template` (`skill://create-template`).

[Switch to 中文](README.zh-CN.md)
