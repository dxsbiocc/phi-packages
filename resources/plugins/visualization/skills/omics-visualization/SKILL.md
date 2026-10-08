---
name: omics-visualization
description: >-
  Render, adapt, and audit publication-ready omics figures from result tables and matrices with bundled templates. Nature-style: claim, hierarchy, restraint, final-size readability. Choose, create, revise, or export scatter, heatmap, bar, distribution, line, network, flow, hierarchy, clustering dendrogram, pathway-ring, gene-structure, or ideogram for transcriptomics, proteomics, metabolomics, enrichment, survival, network, or related omics. Trigger: 组学绘图、科研配图、论文图表、聚类树、通路环图、应力布局,SVG散点,图标柱状图,多通路富集热图,通路分组富集热图,富集点图,GO点图,表达热图,免疫相关热图,关联点阵,平行集合图,泰森多边形,条形甜甜圈,环状比例图,免疫棒棒糖,环形棒棒糖,基因结构,核型图,感兴趣基因,火山图放大,UpSet,韦恩图,Venn,饱和突变热图,ΔΔG热图,oncoprint,突变瀑布图,基因组变异热图,分块聚类热图,NMF相关热图,环状热图,相关气泡热图,对角线分割热图,环形UMAP,UMAP环图,三元图,三元相图,相关散点矩阵,ggpairs,pairs plot,泳道图,游泳图,脊线图,山脊图,joyplot,RidgePlot,基因组覆盖度,coverage track,locus browser,HiChIP,突变棒棒糖,蛋白棒棒糖,g3viz,MutationMapper,蒲公英图,dandelion,Circos,基因组环图,共线性,synteny,弦图,chord diagram,基因组热图,嵌套缩放,nested circos. Not for upstream analysis, business dashboards, AI-generated illustrations, or template libraries.
phi:
  environment: phi:r@1
  attachTo:
  - Visualization
  scripts:
  - name: examples
    description: Find real preview.png images already shipped with omics-visualization templates. Use when the user asks to see examples or styles without providing a data table. This reads the installed catalog and images; it does not create sample data, simulate a plot, render a new figure, or write to the project. Embed the returned preview_markdown images unchanged and explain that they are template examples, not plots of the user's data.
    run:
    - python
    - ./scripts/viz.py
    - examples
    args:
      type: object
      additionalProperties: false
      required:
      - purpose
      properties:
        purpose:
          type: string
          description: Chart family or visual purpose to preview.
        top:
          type: integer
          minimum: 1
          maximum: 4
          default: 4
    approval: read
    output: ./schemas/viz-examples.json
    timeoutSeconds: 120
  - name: route
    description: For new figure selection, profile a result table and shortlist bundled omics templates that fit its purpose. Returns columns and up to six candidates with template_id, fit, risks, and preview_markdown. Use it before choosing a new template; never invent template ids. Do not call it for a small revision of an existing prepared script. To show previews, embed each candidate's preview_markdown and stop for the user's choice. Data and scripts must be inside the project; copy external data in first.
    run:
    - python
    - ./scripts/viz.py
    - route
    args:
      type: object
      additionalProperties: false
      required:
      - data_path
      - purpose
      properties:
        data_path:
          type: string
          format: input-path
          description: Path of the CSV/TSV result table inside the project (project-relative or absolute). It is only read.
        purpose:
          type: string
          description: What the figure should show, for example "volcano plot of differential expression".
        mode:
          type: string
          enum:
          - preview
          - publication
          default: preview
          description: '"publication" is stricter and expects the catalog entry to be confirmed by hand.'
        top:
          type: integer
          minimum: 1
          maximum: 6
          default: 4
        sidecar_dir:
          type: string
          format: input-path
          description: Directory holding companion tables such as nodes.tsv or links.tsv, when not beside the data.
    approval: read
    output: ./schemas/viz-route.json
    timeoutSeconds: 120
  - name: prepare
    description: For a newly selected template, copy its plot.R into the project and return the input contract and editable CONFIG/DATA PREPARATION sections. Edit the copy, then call viz_render. For a revision of an existing prepared plot.R, read and edit that project copy directly; do not copy or reset the bundled template. An existing copy is kept with its edits unless reset is true. Data and scripts must be inside the project; copy external data in first.
    run:
    - python
    - ./scripts/viz.py
    - prepare
    args:
      type: object
      additionalProperties: false
      required:
      - template_id
      - workdir
      properties:
        template_id:
          type: string
          description: A template_id from viz_route, for example "scatter-volcano".
        workdir:
          type: string
          format: project-path
          description: Directory for the copy, inside the project, for example "visualizations/volcano".
        reset:
          type: boolean
          description: Discard an earlier copy and its edits and start again.
    approval: write
    output: ./schemas/viz-prepare.json
    timeoutSeconds: 120
  - name: render
    description: 'Run a prepared plot.R on the input table(s) and write the figure, then check the file. The script, the inputs, and the output must be inside the project; copy external data in first. Returns the output path, format, size in pixels, and the QA result (failed check names). A missing R package or an R error comes back as a short message; do not install packages. A passing QA only means the file is sound: still look at the figure at final size for clipped labels, legends and misleading encodings.'
    run:
    - python
    - ./scripts/viz.py
    - render
    args:
      type: object
      additionalProperties: false
      required:
      - script
      - inputs
      - output
      properties:
        script:
          type: string
          format: input-path
          description: Path of the prepared plot.R (from viz_prepare).
        inputs:
          type: array
          items:
            type: string
            format: input-path
          description: Paths of the table(s) inside the project that the template takes, in the order viz_prepare listed them.
        output:
          type: string
          format: project-path
          description: 'Figure file to write: .png, .pdf or .svg, inside the project.'
        timeout_seconds:
          type: integer
          minimum: 1
          maximum: 900
          default: 180
    approval: write
    output: ./schemas/viz-render.json
    timeoutSeconds: 900
---

# Omics Visualization

Use bundled templates for existing omics results. Make a pattern, comparison,
or uncertainty easier to inspect; select from catalogs, not memory. The chosen
template determines the implementation language.

Default to a Nature-style figure: one defensible claim, restrained color and
readable type at final size, unless journal rules specify otherwise.

## Installed examples without data

Without user data, call `viz_examples` by purpose. Embed its existing
`preview.png` via `preview_markdown` unchanged as a Markdown image; label it a
template example. Show up to four; ask for choice or data. Do not create
simulated data, render, copy, or replace it; report no match.

## Revise an existing figure

For a prior figure, follow [revision workflow](references/revision-workflow.md):
edit its script and palette; do not call `viz_route` or `viz_prepare` for colors.

## Match a reference image

For a new figure guided by a user image, follow
[reference workflow](references/reference-figure-workflow.md). Inspect the image
before routing. If editable source already exists, revise it instead.

## New figure route

Use `viz_route`, `viz_prepare`, `viz_render` for new figures.

1. State the visualization purpose before choosing a chart: what should become
   easier to see, compare, verify, or question after plotting? If the purpose is
   unclear, infer the smallest honest purpose from the user's request and data;
   ask only when several purposes would require materially different figures.
2. Inspect the data's observation unit, columns, comparison, and claim. Check
   row count, types, category order, ranges, missingness, zeros, duplicate IDs,
   and any pairing, time, hierarchy, interval, network, or matrix structure.
   `viz_route` gives column types; inspect values when the figure depends on them.
3. Write or infer a one-sentence figure claim before selecting geometry. A
   single composite glyph (for example a circular tree inset in a polar track)
   is still one template. Follow **Multi-panel composition** only for labelled
   a/b/c figures that combine several plots with distinct evidence roles.
4. Shortlist templates with `viz_route` (data path, purpose, mode).
   It profiles the table against
   [references/template_contracts.json](references/template_contracts.json) and
   scans recognized companion files in the input directory (`nodes.tsv`,
   `links.tsv`, `rowInfo.tsv`, `colInfo.tsv`, `enrichment.tsv`, `cytoband.tsv`,
   `domains.tsv`, `karyotype.tsv`; pass `sidecar_dir` when they live elsewhere),
   reporting whether their IDs align. Alignment warnings lower confidence; they
   do not authorize silent filtering, reordering, imputation, or upstream
   statistical analysis. A high-confidence result is a shortlist, not permission
   to skip step 6.
5. If the router is low confidence, if several viable templates would change
   the scientific reading, or if the output is publication-critical, read
   [references/plots.yaml](references/plots.yaml) and then only the relevant
   family catalog under `references/catalog/`. Compare `use_when`, `avoid_when`,
   and `input_shape`. When two templates of one family look alike, read
   [references/mis-routes.md](references/mis-routes.md). If one template clearly
   matches the purpose, distribution, claim, and data contract, proceed with it
   and record the reason. If alternatives change interpretation, show up to
   four real previews (`preview_markdown`), fit, risks and required input;
   let the user choose.
6. Select one canonical template and make a project-local copy with
   `viz_prepare`. It returns the template's purpose, the tables it takes, its R
   dependencies, and the CONFIG and DATA PREPARATION sections with their line
   numbers; read the PLOT lines only when the geometry must change. [Without the
   tool: read the template's `plot.R` in full, and any helper it sources, before
   modifying or executing it. Current templates use
   [scripts/lib/common.R](scripts/lib/common.R).]
7. Pick colors from the generated [palette index](references/palettes.yaml).
   Copy hex from that file or from
   [full color catalog](references/palettes/colors.json). These are byte-for-byte
   copies of Phi's shared `resources/palettes/` catalog so the plugin remains
   self-contained; edit only the shared source when maintaining the plugin.
   Do not invent hex codes. Prefer `Qualitative.Safe` for discrete groups,
   `Quantitative.BluGrn` for heatmaps, and `Diverging.RdBu` for signed values,
   unless the user names a palette or asks to match OmicsAgent (`Brand.Algolia`)
   or an existing figure. Reserve saturated color for the focal comparison and
   keep secondary marks neutral when the template permits it.
8. Keep the working copy project-local: under the
   active Phi project working directory, not the input data file's parent
   directory. Treat data directories outside the
   project as read-only inputs even if the shell can technically write there.
   Keep generated figures, adapted source files, scratch files, reports, and QA
   artifacts under the project working directory so Phi can preview/open them.
   If no output directory is specified, create a concise directory such as
   `visualizations/<short-task-name>/` or `plots/<short-task-name>/` inside the
   project. `viz_prepare` and `viz_render` refuse paths outside the project, and
   `viz_prepare` points the copy at the installed read-only `scripts/lib/common.R`
   by absolute path and returns that location as `common_r`. Preserve the
   prepared bootstrap when adapting a plot; use `common_r` for any helper
   reference instead of guessing repository resource paths. `viz_prepare` and
   `viz_render` refresh missing installed helper references in existing copies
   without resetting their plotting code. [By hand: source that helper by absolute path, or copy only
   the helper and set `OMICS_VISUALIZATION_SKILL_ROOT` to the installed skill
   root]. Do not copy `references/`, `references/palettes/`, or catalog files into
   the project merely to make palette lookup work; choose palette ids from the
   installed generated catalog and leave it in the skill. Do not
   edit the installed skill copy during an ordinary plotting task.
9. Adapt the visible source directly:
   - edit `CONFIG` for column mappings, labels, and ordinary presentation
     (`circular`, `outer = "bar"` vs `"point"`, hole, inset, SVG file).
     A layout flag or a different glyph file is not a new template id.
   - edit `DATA PREPARATION` when the input shape differs;
   - edit `PLOT` for structural visual changes.
   Do not introduce an external plotting configuration DSL.
   Do not add a sibling template that only changes `layout`, `circular`,
   an algorithm alias (kk / fr / lgl), or the SVG glyph file.
10. Set the final output size before polishing. For manuscript-like output,
   read [references/nature-figure-principles.md](references/nature-figure-principles.md)
   and apply only the parts relevant to the selected template.
11. Render with `viz_render` (the prepared script, the input table(s), an output
   `.png`, `.pdf` or `.svg`) [run `Rscript plot.R <input> <output>` as documented
   at the top of the script]. If dependencies are missing, report them; do not
   install packages without authorization and do not silently switch
   implementations.
12. `viz_render` runs the lightweight artifact QA and names the checks that
   failed. That
   only shows the file is sound: inspect the rendered artifact at final size and
   correct labels, scales, legends, clipping, overlaps, spacing, and misleading
   encodings, then rerun until the output and source agree.
Routing is complete only when the chosen template, its input assumptions,
implementation language, sourced helpers, and runnable project-local copy are
known. When the contract router is used, delivery notes should include the
recommended template id, confidence, decisive matched shape, and any listed
risk that had to be checked. Delivery is complete only after the final rendered
artifact has passed lightweight QA and visual inspection.

## Efficient use

- Use `viz_route` as the fast path for preview-mode and
  common result-table shapes. Trust high-confidence contract recommendations
  enough to avoid reading unrelated family catalogs, but still confirm the
  selected template's contract (`viz_prepare`) before rendering.
- Trust `references/plots.yaml` and the family catalogs for low-confidence,
  publication, or novel shapes. Do not search the template scripts by text
  unless a catalog entry is missing, stale, or internally inconsistent.
- When adding or editing routing contracts, run
  `skill_run({ skill: "omics-visualization", script: "validate_template_contracts.py", args: ["--json"] })` before delivery.
  The validator must pass with no errors; warnings require an explicit note.
  Use the text report with `--max-missing-per-family 3` to prioritize future
  contract coverage without manually scanning the whole catalog.
- When the user names an exact template id, or the data shape leaves one clear
  candidate, proceed without a candidate-selection pause.
- Do not render every candidate. Use shipped previews for material alternatives,
  then run only the chosen project-local source.
- Prefer `CONFIG` edits, label tightening, scale changes, and legend pruning
  before touching the `PLOT` section. Structural rewrites are for misleading
  geometry, not ordinary polish.
- Keep final QA proportional: single plots get `viz_render`'s QA plus a
  final-size visual inspection; labelled multi-panel figures also get the
  Patchwork geometry audit.

## Nature-style figure defaults

- Argument first: every figure should answer one claim or one bounded
  comparison. Remove panels, labels, encodings, and legends that do not help a
  reader evaluate that claim.
- Restraint over decoration: no ornamental gradients, busy backgrounds,
  gratuitous icons, or palette changes that encode nothing. Use whitespace and
  alignment as the main organizing devices.
- Typography is judged at final size. Use concise sentence-case labels with
  units; avoid local panel titles when the caption can carry the narrative.
  As a practical default, keep tick labels at least 6 pt, axis and legend text
  at least 7 pt, and panel letters 8-10 pt bold.
- Hierarchy should be visible before details. Give the decisive evidence the
  most space or clearest contrast; make controls, references, and secondary
  strata quieter but still legible.
- Encodings must be biologically and statistically honest. Prefer effect size
  plus uncertainty when supplied; avoid summary bars for distributions when a
  raw-point, interval, boxplot, violin, raincloud, or ridge template better
  exposes the data shape.
- Template choice follows the data. Do not force a requested chart type when
  the distribution snapshot shows it would hide sample size, censoring,
  sparsity, outliers, paired structure, or uncertainty that the figure purpose
  depends on; recommend a better-fitting catalog template and explain the
  tradeoff.
- Color must survive reduction and color-vision checks. Use catalog palettes,
  avoid rainbow/jet and pure red-green dependence, and add redundant shape,
  linetype, label, or ordering when color carries a critical distinction.

## Confusable templates

Read [references/mis-routes.md](references/mis-routes.md) when two templates of one
family look alike or the router's pick and the user's wording disagree. It lists
the confusable pairs and the geometry that separates them. Confirm against the
family catalog either way.

## Multi-panel figures

For a labelled multi-panel figure, figure assembly, or rearrangement of several
plots, follow [references/multipanel-workflow.md](references/multipanel-workflow.md).
Its delivery bar: every panel has a distinct evidence role, the rendered geometry
audit passes, and the final-size visual inspection finds no unresolved clipping,
collision, hierarchy, or reading-order defect.

## Language boundary

Do not search for or create a duplicate bundled implementation in another
language. If the user explicitly requires a language not represented by the
selected canonical template, state that limitation and create a project-local
implementation only when it is necessary to satisfy the request; do not add it
to the bundled library as a second canonical template.

## Scientific integrity

- For volcano plots, state the adjusted-p-value and absolute log2-fold-change cutoffs. Use the same values for point classification, cutoff lines, legend text, and any reported Up/Down counts. Keep padj-only differential-expression counts distinct from counts that also apply a fold-change cutoff; neither threshold is universal, so honor the user's values or disclose the selected template's values.
- Do not invent sample sizes, statistical tests, p-values, adjusted p-values,
  effect sizes, uncertainty, group mappings, or biological interpretations.
- Do not silently filter, aggregate, impute, clip, coerce, or sample data.
  Report every material transformation and its before-and-after row count.
- Preserve identifiers and biological units. Flag duplicated identifiers,
  missing values, non-finite values, invalid ranges, and ambiguous columns
  before plotting.
- Distinguish raw from adjusted p-values and technical from biological
  replicates. Do not infer biological importance from visual or statistical
  separation alone.
- If the figure requires upstream analysis that has not been performed,
  identify that prerequisite instead of fabricating results.
- Display clustering (dist + hclust on a supplied matrix) is allowed when
  that is the figure (`tree-dendrogram`, clustered heatmaps). Do not present
  it as an independent statistical result, and do not invent the matrix.
- Use catalog palettes. Do not invent hex. Do not interpolate hex when a
  palette is too short; `palette_colors()` appends unused Qualitative then
  Brand colours. Do not use Artwork or Concept palettes unless the user
  asks for a decorative theme.
- Plot canvases default transparent; dark viewers may show black alpha. For a
  newly rendered figure with hidden dark labels, export a white-backed PNG
  preview and retain vector PDF/SVG. Never re-render shipped example PNGs.

## Delivery

Return the final rendered figure and the exact project-local source used to
create it. Include the visualization purpose, selected template ID, selection
rationale, input-to-column mapping, material transformations, and unresolved
limitations in a concise handoff. Prefer SVG or PDF for editable scientific
graphics and add a PNG preview when useful. Do not claim rendering or visual
validation when either step was blocked.

## Scope

This skill applies existing templates to real data and revises or audits the
resulting figures. Creating, cataloging, testing, or maintaining reusable
templates belongs to a template-authoring workflow, not this skill.
