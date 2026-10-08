# Easy mis-routes

Read this when two templates of one family look alike, or when the router's
recommendation and the user's wording disagree. Paths and template ids are the
same as in SKILL.md.

Confirm against the family catalog. These pairs share vocabulary but not
geometry:

- Ranked enrichment list → `bar-enrichment-*`. Classic GO/KEGG
  dotplot (x = gene ratio, size = Count, fill = −log10(p)) →
  `bar-enrichment-dot`. Facet by ontology is CONFIG. Supplied class
  → term tree with Count and p on a circumferential track →
  `tree-enrichment-ring`.
  That figure is one template, not a multi-panel layout. Polar track plus
  circular dendrogram are combined with inset_element; do not treat them
  as two catalog ids. A DE expression heatmap whose right-hand panels are
  enrichment bars aligned to Up/Down gene splits → `heatmap-enrichment-zoom`.
  Extra databases are extra rows in the enrichment table, not extra ids.
  The same terms across contrasts, fill = signed −log10(p), missing
  cells marked not tested → `heatmap-enrichment-terms`. Do not run
  enricher or convert species in the script. That is not
  heatmap-enrichment-zoom and not heatmap-corr-dot.
  Rectangular association of two category axes with supplied rho and p,
  colour = rho, shape = p cutoff → `heatmap-corr-dot`. Do not compute
  Spearman in the script. That is not heatmap-signif (square matrix
  from samples) and not heatmap-two (two variable sets from samples).
  Colour = rho, size = −log10(p) or a supplied magnitude, stars from
  supplied p → `heatmap-corr-bubble`. Facet columns or vline
  intercepts stay in CONFIG. That is not heatmap-corr-dot (shape
  encodes the p cutoff). Do not compute Spearman in the script.
  Each cell split on the diagonal, two continuous fills (supplied
  coefficient and p) → `heatmap-corr-triangle`. Optional p-dots
  stay in CONFIG. That is not heatmap-corr-dot or
  heatmap-corr-bubble. Do not compute Spearman or impute missing
  pairs in the script.
  Site × amino-acid ΔΔG tiles with a dual-y mean overlay →
  `heatmap-mutation-energy`. That is one template, not heatmap-basic
  plus a line chart. Do not interpolate tiles or estimate energies
  in the script. Mean ΔΔG is the mean of the supplied cells at that
  site. A protein-residue lollipop with class pies is
  `ideogram-lollipop`, not this heatmap. Gene × sample layered
  alteration glyphs (CNV fill, mutation bar, frequency bars) →
  `heatmap-oncoprint`. That is not heatmap-mutation-energy, not
  ideogram-lollipop, and not bar-waterfall. Glyph and bar side
  stay in CONFIG. Do not parse MAF or fetch cBioPortal in the
  script.
  A supplied rectangular score or Pearson matrix clustered with
  dendrograms cut into k blocks and an optional row-group strip →
  `heatmap-cluster-block`. That is not heatmap-cluster-basic (those
  rows are z-scored samples). Do not compute correlation in the
  script. cutree k stays in CONFIG.
  A circular heatmap whose sectors are a supplied Group, with optional
  q-value diamonds and an Euler hole → `heatmap-circos-split`. That is
  one template, not heatmap-basic plus a Venn. inner and the q-value
  track are CONFIG. Do not cluster or rescale in the script.
- Named parent/path hierarchy → `tree-basic`. Numeric matrix clustered with
  dist + hclust → `tree-dendrogram`. `circular = TRUE` stays in CONFIG.
- Supplied node x/y → `graph-force`, `graph-basic`, or `graph-cartesian`.
  Edge table with a computed layout → `graph-stress`. `graph-force` uses
  `layout = "manual"`; it does not run a force algorithm.
  Equal-spacing single-ring nodes → `graph-circular`. Hub nodes on an
  inner ring plus category-grouped outer nodes with class-coloured
  weighted edges → `graph-circular-concentric`. That is not
  graph-circular and not graph-chord.
- Sunburst when ring **angle** encodes descendant size. Outer bar height or
  point size is not a sunburst. A circos stacked doughnut of within-group
  composition with radial bars on the inner track → `pie-doughnut-bar`.
  That is one template; sector width and item labels stay in CONFIG.
- Ordinary pch scatter → `scatter-group`. SVG glyphs at x-y → `scatter-svg`.
  A supplied UMAP in a circos whose sectors are cell types, with
  stacked metadata rings → `scatter-umap-circos`. Extra tracks and
  `scale = "log10"` vs `"linear"` stay in CONFIG. Do not run UMAP or
  wrap plot1cell. A doughnut without embedding is `pie-doughnut-bar`.
  Three non-negative components as a composition (one point in a
  simplex) → `scatter-ternary`. That is not scatter-group. Discrete vs
  continuous colour, and optional size, stay in CONFIG. Do not compute
  lineage scores, NMF, relative abundance, enrichment, or absorption
  probabilities in the script. Closing a + b + c to 1 is a display
  transform. Ternary KDE / percentile density bands are not this id.
  Several numeric columns as a pairs grid (lower scatter + lm,
  diagonal names, upper Pearson r tiles) → `scatter-pairs`. That is
  not scatter-matrix (category × category bubbles) and not
  scatter-correlation (one pair). More than about 12 variables →
  `heatmap-signif`. Method and gap stay in CONFIG. Do not run DE or
  enrichment in the script.
  Category SVG instead of bar y-axis text → `bar-svg-icon`. Changing the
  SVG file is CONFIG, not a new id. A signed delta lollipop (stem through
  zero, point size = -log10(p), fill = signature set) →
  `scatter-lollipop-delta`. That is not scatter-cleveland and not a
  dumbbell. Do not add an id per immune deconvolution method.
  Grouped circular lollipops with sector fans, radial error bars, and
  outer group labels → `scatter-lollipop-circular`. That is not
  scatter-lollipop-polar (plain polar stems) and not
  scatter-lollipop-radial (a fan). se is supplied; do not compute it
  in the script. Subgroup names sit in the points. Point radius is
  the mean; whiskers are se. Do not print numeric mean labels.
  Mutations on a protein amino-acid axis (stems, pie or circle heads,
  domain bar) → `ideogram-lollipop`. That is not scatter-lollipop-delta,
  scatter-rank, or scatter-cleveland. Pie vs circle is CONFIG. Do not
  parse MAF, fetch cBioPortal, or look up Pfam in the script. A gene ×
  sample oncoprint is `heatmap-oncoprint`, not ideogram-lollipop.
  Dense sites collapsed into clustered dandelions (stem = cluster
  size or mean score; fan / pie / circle / pin heads) →
  `ideogram-dandelion`. That is not ideogram-lollipop (one stem per
  site) and not ideogram-density. type and maxgaps are CONFIG. Do
  not fetch TxDb, UCSC, or VCF in the script.
- Transcript exon / CDS / UTR tracks → `ideogram-gene`. Chromosome
  G-banding with no loci → `ideogram-karyotype`. Window statistic fill
  along a chromosome → `ideogram-density`. A supplied gene-loci table
  (points + names, colour by category) → `ideogram-loci`. Catalog
  preview is circular. Vertical vs circular is `config$layout`, not a
  second id. `coord_flip` and hg19 vs hg38 stay in CONFIG / the
  cytoband sidecar. Stacked coverage on one genomic window (scATAC /
  HiChIP / ChIP signal, optional loop arcs and gene bars) →
  `ideogram-coverage`. Extra tracks, a highlight, and loops on/off
  are CONFIG / sidecars. That is not ideogram-density, not
  ideogram-gene, and not line-*. Do not call peaks, loops, MACS, or
  Hi-C in the script. Protein mutation lollipops are
  `ideogram-lollipop`, not ideogram-gene. Dense clustered dandelions
  are `ideogram-dandelion`, not ideogram-lollipop.
  One genome as circular sectors with concentric tracks
  (histogram / scatter / line) and optional SV links →
  `ideogram-circos`. Extra tracks and links on/off are CONFIG.
  Do not call CNV or density in the script. Two assemblies with
  homology ribbons → `ideogram-synteny`. That is not
  ideogram-circos and not graph-chord. Do not run BLAST in
  the script. A Circos chord whose sector width is total flow
  → `graph-chord`. Equal-spacing circular nodes are
  graph-circular. Genomic links on an ideogram are
  ideogram-circos. Genomic interval tiles with connector
  lines → `ideogram-circos-heatmap`. That is not
  heatmap-circos-split. side inside/outside is CONFIG. Do
  not call peaks or z-score in the script. An outer genome
  plus inner zoom windows joined by correspondence →
  `ideogram-nested`. Extra windows are extra rows. Do not
  call DMR in the script. That is not ideogram-circos.
- Signed network: colour encodes sign of supplied corr, width encodes
  `|corr|`. Do not compute correlation in the script. Label size is fixed
  and readable; do not bind it to `|corr|` or tile area.
  Several categorical columns plus a supplied count as parallel-axis
  ribbons → `sankey-parallel-sets`. That is not sankey-basic or
  sankey-level (those take a source-target edge list). Do not tabulate
  raw samples in the script.
  Supplied x-y points filled as nearest-site tiles → `scatter-voronoi`.
  That is not scatter-group (points only) and not scatter-contour
  (convex hull). Do not run PCA, UMAP, or clustering in the script.
  Tile area is not density.
  Overview volcano → `scatter-volcano`. Rectangular zoom plus a densely
  labelled inset of that window → `scatter-volcano-inset`. That is one
  template, not a multi-panel a/b layout and not a CONFIG flag on
  scatter-volcano. Do not run DE in the script.
  Combination membership with top intersection-size bars and left set-size
  bars → `bar-upset`. Two bar series (observed vs expected) and matrix
  colour by n-sets stay on this id. That is one template, not a
  multi-panel layout. Do not tabulate raw samples or compute expected
  counts in the script. Overlapping circles of 2–4 supplied sets →
  `bar-venn`. Percentage and element-name labels stay in CONFIG. Five
  or more sets stay on `bar-upset`. Do not add an id per ggvenn / venn
  package or for a 5–7 petal Venn.
  One lane per patient on a time axis, treatment intervals as
  rectangles, response points and stop marks → `bar-swimmer`. Extra
  treatments and left-hand annotation tiles stay in CONFIG. That is
  not bar-stack, not scatter-dumbbell, and not line-survival. Do not
  convert dates or call response in the script.
  Stacked 1D density ridges by group, optionally faceted by feature →
  `line-ridge`. That is not boxplot-violin (mirrored density) and not
  a raincloud. Marker facets and left group bars stay in CONFIG. Do
  not cluster cells or pool an All-cells row in the script; that
  group is supplied. The vertical tick is a display quantile of the
  supplied values. Genome-coordinate coverage tracks are
  `ideogram-coverage`, not line-ridge.
  Normal vs Tumor (or any two-level contrast) across many cancer
  types with alternating background bands and per-category stars →
  `boxplot-differential-bg`. That is not boxplot-group and not
  boxplot-differential-expression. show_ns stays in CONFIG. Do not
  filter matched patients or run DE in the script.

