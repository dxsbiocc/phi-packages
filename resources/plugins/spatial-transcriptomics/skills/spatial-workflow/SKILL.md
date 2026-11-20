---
name: spatial-workflow
description: >-
  Run a spatial transcriptomics workflow on 10x Visium (incl. HD bin/seg) or BGI
  Stereo-seq (GEM/GEF, square bin or cellbin): inspect platform/bin, QC,
  normalization, variable genes, PCA, Leiden clustering, UMAP, the
  tissue-coordinate spatial graph, spatial cluster/gene maps, marker genes,
  panel cell-type annotation, spatially variable genes (Moran's I / Geary's C),
  neighbourhood enrichment and co-occurrence, deconvolution against a
  single-cell reference (scipy NNLS), ligand-receptor communication (squidpy
  ligrec), and cross-sample differential spatial patterns. Use for 空间转录组,
  spatial transcriptomics, Visium, Visium HD, Stereo-seq, 华大, squidpy,
  spatialdata, 空间可变基因, SVG, 邻域富集, deconvolution, 细胞通讯, .gef, .gem,
  or when a spatial matrix needs a reproducible analysis chain. Python route
  (scanpy 1.12 + squidpy 1.8 + spatialdata); the R/OSTA route is a future
  sibling. Not for raw-read processing (Space Ranger / SAW are wrappers) or
  scATAC-seq.
phi:
  environment: plugin:spatial
  attachTo:
    - SpatialTranscriptomics
  scripts:
    - name: py_inspect
      description: >-
        Read a Visium Space Ranger directory or .h5, a Stereo-seq GEF/GEM
        (square bin or cellbin), a Visium HD bin table (.csv/.tsv), or an .h5ad
        with spatial coordinates, and report its platform, bin, shape, layers,
        embeddings, whether spatial coordinates are present, and the metadata
        columns worth using. Lists which follow-up steps the object is ready
        for. Call this first on an unfamiliar input. With output, saves the
        loaded object as .h5ad so later steps can chain.
      run: [python, ./scripts/project_inspect.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A Visium directory/.h5, a Stereo-seq .gef/.gem, a Visium HD bin table, or an .h5ad with spatial coordinates.
          output:
            type: string
            format: project-path
            description: Where to save the loaded AnnData object as .h5ad, inside the project.
          platform:
            type: string
            enum: [auto, visium, visium-hd, stereoseq]
            description: Force a platform; auto detects from the input (default auto).
          bin_size:
            type: integer
            description: Stereo-seq square bin / Visium HD bin side; null reads the bin the file carries. Choosing a bin is a biological decision.
          cellbin:
            type: boolean
            description: Read Stereo-seq cellbin (segmented cells) instead of square bins; segmentation is done upstream by SAW.
          min_units:
            type: integer
            minimum: 0
            description: Drop spatial units with fewer genes.
          min_genes:
            type: integer
            minimum: 0
            description: Drop genes detected in fewer units.
          var_names:
            type: string
            enum: [gene_symbols, gene_ids]
            description: Which 10x feature column becomes var_names (default gene_symbols).
      approval: read
      timeoutSeconds: 600
    - name: py_qc
      description: >-
        Compute mitochondrial, ribosomal and hemoglobin fractions per spatial
        unit, apply the given thresholds, and report how many units each
        threshold removed plus the metric distributions behind a sensible
        choice. Writes before/after violin panels and the same metrics mapped
        onto tissue. Spatial units are spots/bins/cells, not droplets, so there
        is no doublet gate. Start from raw counts.
      run: [python, ./scripts/qc_filter.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A spatial input (Visium/Stereo-seq/.h5ad) with raw counts.
          output:
            type: string
            format: project-path
            description: Where to save the filtered object as .h5ad.
          figures_dir:
            type: string
            format: project-path
            description: Directory for QC plots; omit to skip figures.
          platform:
            type: string
            enum: [auto, visium, visium-hd, stereoseq]
            description: Force a platform (default auto).
          bin_size:
            type: integer
            description: Stereo-seq / Visium HD bin side; null reads the file's bin.
          cellbin:
            type: boolean
            description: Read Stereo-seq cellbin instead of square bins.
          min_genes:
            type: number
            description: Minimum genes per spatial unit (default 100; lower than single-cell because units are shallower).
          min_counts:
            type: number
            description: Minimum UMIs per unit; below this a unit is likely empty.
          max_counts:
            type: number
            description: Maximum UMIs per unit; null disables the upper gate.
          mt_max:
            type: number
            description: Maximum mitochondrial percent (default 20; spatial tissue tolerates more than dissociated cells).
          mt_pattern:
            type: string
            description: Regex for mitochondrial genes (default ^MT- for human; use ^mt- for mouse).
          quantiles:
            type: array
            items:
              type: number
            description: Quantiles to report per metric (default 0.01 and 0.99).
          sample_column:
            type: string
            description: obs column naming each unit's sample, for a per-sample breakdown.
      approval: write
      timeoutSeconds: 1800
    - name: py_preprocess
      description: >-
        Normalize (log1p of library-size-normalized counts), select variable
        genes on the counts layer (seurat_v3), snapshot expression to lognorm,
        scale, run PCA, build the expression neighbor graph, cluster with Leiden
        at one or more resolutions, embed with UMAP, then build the spatial
        adjacency graph from tissue coordinates and map clusters onto tissue.
        Reports variance explained, cluster sizes, top variable genes and the
        spatial graph density. Reads an already-QC'd object.
      run: [python, ./scripts/preprocess_cluster.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: AnnData saved as .h5ad, after st_py_qc.
          output:
            type: string
            format: project-path
            description: Where to save the processed object as .h5ad.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the UMAP, elbow, HVG and spatial cluster plots.
          scale_factor:
            type: number
            description: Library-size target sum (default 10000).
          n_hvg:
            type: integer
            description: Number of variable genes (default 2000).
          hvg_flavor:
            type: string
            description: HVG flavor on the counts layer (default seurat_v3).
          n_pcs:
            type: integer
            description: Principal components to compute (default 30).
          dims:
            type: integer
            description: PCs used for neighbors and UMAP (default 30).
          resolution:
            type: array
            items:
              type: number
            description: One or more Leiden resolutions, for example 0.8 (repeat the flag for several).
          run_umap:
            type: boolean
            description: Run UMAP (default true).
          neighbor_method:
            type: string
            enum: [knn, radius, delaunay, grid]
            description: Spatial adjacency geometry (default knn); Visium hex suits grid/delaunay, Stereo-seq bins suit knn/radius.
          n_neighs:
            type: integer
            description: Neighbours for knn/grid adjacency (default 6, the Visium hex layout).
          radius:
            type: number
            description: Radius for radius adjacency in coordinate units; null derives it from spacing.
          coord_key:
            type: string
            description: obsm key holding spatial coordinates (default spatial).
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      timeoutSeconds: 3600
    - name: py_spatial
      description: >-
        Spatial statistics on a clustered, normalized object: Moran's I or
        Geary's C spatially variable genes, neighbourhood enrichment and
        co-occurrence of domains, and (when group_column names a condition)
        cross-sample differential spatial patterns. Rebuilds the spatial
        adjacency graph from tissue coordinates. Reports the SVG ranking, the
        top enriched domain pairs, and genes whose spatial structure changes
        between conditions. A section that needs a missing input reports the
        gap and continues.
      run: [python, ./scripts/spatial_statistics.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: AnnData saved as .h5ad, after st_py_preprocess (clusters + lognorm).
          output:
            type: string
            format: project-path
            description: Where to save the object with spatial statistics attached.
          figures_dir:
            type: string
            format: project-path
            description: Directory for SVG, neighbourhood and differential plots.
          group_by:
            type: string
            description: obs column defining the domains (default leiden).
          coord_key:
            type: string
            description: obsm key holding spatial coordinates (default spatial).
          neighbor_method:
            type: string
            enum: [knn, radius, delaunay, grid]
            description: Spatial adjacency geometry (default knn).
          n_neighs:
            type: integer
            description: Neighbours for knn adjacency (default 6).
          radius:
            type: number
            description: Radius for radius adjacency; null derives it from spacing.
          svg_method:
            type: string
            enum: [moran, geary]
            description: Spatially variable gene statistic (default moran).
          svg_corr_method:
            type: string
            description: Multiple-testing correction, a statsmodels name (default fdr_bh; benjamini-hochberg is an alias).
          n_perms:
            type: integer
            description: Permutations for neighbourhood enrichment (default 1000).
          cooc_interval:
            type: integer
            description: Distance intervals for co-occurrence (default 50).
          group_column:
            type: string
            description: Condition/sample obs column for cross-sample differential; null skips it.
          min_units:
            type: integer
            description: Minimum units for a condition to enter the differential contrast (default 10).
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      timeoutSeconds: 3600
    - name: py_markers
      description: >-
        Rank marker genes per spatial domain with rank_genes_groups on the
        pre-scale lognorm snapshot, write one CSV per domain, and draw a rank
        plot and heatmap. Reports the top genes and log2 fold changes per
        domain. The p-values are descriptive (clusters were defined from the
        same data); treat the ranking as evidence, not significance.
      run: [python, ./scripts/find_markers.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: AnnData saved as .h5ad, after clustering.
          output:
            type: string
            format: project-path
            description: Where to save the object with markers attached.
          tables_dir:
            type: string
            format: project-path
            description: Directory for one markers-<domain>.csv per domain.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the marker rank plot and heatmap.
          group_by:
            type: string
            description: obs column to group by (default leiden).
          test_use:
            type: string
            enum: [wilcoxon, t-test, t-test_overestim_var, logreg]
            description: Ranking test (default wilcoxon).
          logfc_threshold:
            type: number
            description: Minimum log2 fold change (default 0.25).
          min_pct:
            type: number
            description: Minimum fraction of units expressing the gene (default 0.1).
          only_positive:
            type: boolean
            description: Report positive markers only (default true).
          top_n:
            type: integer
            description: Top genes to report per domain (default 10).
      approval: write
      timeoutSeconds: 3600
    - name: py_annotate
      description: >-
        Score every domain against a JSON marker panel (cell type -> gene
        symbols) with a panel-mean-minus-rest rule, and label domains either
        from the highest score or from an explicit mapping JSON (cluster ->
        label). Writes the label column into obs, draws a marker dot plot and a
        spatial cell-type map. A spatial label names the dominant cell type or
        region of a spot/bin, not a single pure cell. Without panel_path or
        mapping_path nothing is annotated (a script must not invent labels).
      run: [python, ./scripts/annotate_clusters.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: AnnData saved as .h5ad, after clustering (and markers for the dot plot).
          output:
            type: string
            format: project-path
            description: Where to save the annotated object as .h5ad.
          panel_path:
            type: string
            format: input-path
            description: JSON object mapping cell type to an array of gene symbols.
          mapping_path:
            type: string
            format: input-path
            description: JSON object mapping every cluster to a label; overrides panel scoring.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the dot plot and spatial cell-type map.
          group_by:
            type: string
            description: Cluster column (default leiden).
          label_column:
            type: string
            description: obs column to write labels into (default cell_type).
          min_score:
            type: number
            description: Minimum panel score to accept a label (default 0).
          coord_key:
            type: string
            description: obsm key holding spatial coordinates (default spatial).
      approval: write
      timeoutSeconds: 3600
    - name: py_deconvolve
      description: >-
        Estimate the cell-type proportions of each spatial unit by non-negative
        least squares (scipy NNLS) against a single-cell reference .h5ad. Writes
        obsm['deconvolution'] and a per-type proportion column, and maps the top
        types onto tissue. This is a linear estimate, NOT a trained model;
        cell2location / tangram-sc / stereopy are pip-only and need a project
        environment. Without reference_path the step reports the gap and
        continues. The reference comes from the single-cell plugin.
      run: [python, ./scripts/deconvolve.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: AnnData saved as .h5ad, with a counts layer and spatial coordinates.
          output:
            type: string
            format: project-path
            description: Where to save the object with proportions attached.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the spatial proportion maps.
          reference_path:
            type: string
            format: input-path
            description: Single-cell reference .h5ad; null skips deconvolution with a message.
          reference_cell_type_key:
            type: string
            description: obs column in the reference naming cell types (default cell_type).
          deconv_method:
            type: string
            enum: [nnls]
            description: Deconvolution solver (nnls in-plugin); model-based tools need a project environment.
          coord_key:
            type: string
            description: obsm key holding spatial coordinates (default spatial).
          neighbor_method:
            type: string
            enum: [knn, radius, delaunay, grid]
            description: Spatial adjacency geometry (default knn).
          n_neighs:
            type: integer
            description: Neighbours for knn adjacency (default 6).
      approval: write
      timeoutSeconds: 3600
    - name: py_communication
      description: >-
        Infer ligand-receptor cell-cell communication between domains with
        squidpy ligrec (a CellPhoneDB-style permutation test). Runs offline with
        a curated starter ligand-receptor set unless interactions_path supplies a
        fuller table (source/target or ligand/receptor columns). Reports
        significant pairs per domain combination. liana / CellPhoneDB core are
        pip-only and need a project environment. Needs clusters and the lognorm
        snapshot.
      run: [python, ./scripts/communication.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: AnnData saved as .h5ad, after st_py_preprocess (clusters + lognorm).
          output:
            type: string
            format: project-path
            description: Where to save the object with communication results attached.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the ligrec dot plot.
          group_by:
            type: string
            description: obs column defining the domains (default leiden).
          interactions_path:
            type: string
            format: input-path
            description: Ligand-receptor TSV/CSV (source/target or ligand/receptor); null uses the shipped offline set.
          ccc_threshold:
            type: number
            description: p-value threshold for significant interactions (default 0.01).
          n_perms:
            type: integer
            description: Permutations for the ligrec test (default 1000).
          coord_key:
            type: string
            description: obsm key holding spatial coordinates (default spatial).
          neighbor_method:
            type: string
            enum: [knn, radius, delaunay, grid]
            description: Spatial adjacency geometry (default knn).
          n_neighs:
            type: integer
            description: Neighbours for knn adjacency (default 6).
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      timeoutSeconds: 3600
    - name: py_notebook
      description: >-
        Render a runnable spatial analysis notebook into the project from the
        fragment library, then open it with Phi's notebook tools and run cells
        on demand. The template is rendered, not copied: placeholders are filled
        from the input data and the requested steps, and every tunable value
        lands in one config cell at the top that the analysis cells reference.
        Pick a preset with module (1 QC, 2 preprocess+cluster, 3 spatial domains
        +markers+annotation, 4 SVG+neighbourhoods, 5 deconvolution+communication,
        6 full chain) or name steps directly; dependencies are added
        automatically. inject adds an external Python script as a code cell;
        merge splices in another notebook, whole or a cell range. Starts from a
        Visium directory/.h5, a Stereo-seq GEF/GEM, a Visium HD bin table or an
        .h5ad: raw-read processing (Space Ranger, SAW) and cell segmentation are
        upstream wrappers, not notebook cells. Use self_check to execute every
        cell in the output directory now and fail loudly if one breaks.
      run: [python, ./scripts/build_notebook.py]
      args:
        type: object
        additionalProperties: false
        required: [output, data]
        properties:
          module:
            type: integer
            enum: [1, 2, 3, 4, 5, 6]
            description: Spatial preset. Pass module or steps, not both.
          steps:
            type: array
            items:
              type: string
              enum: [setup, load, qc-metrics, qc-plots, filter, normalize, hvg, scale, pca, cluster, spatial-neighbors, spatial-plots, markers, annotate, assign-labels, svg, nhood, deconvolve, communication, spatial-diff, save]
            description: >-
              Steps to render instead of a preset; declared dependencies are
              added automatically. Pass module or steps, not both.
          output:
            type: string
            format: project-path
            description: Path of the .ipynb to write inside the project.
          data:
            type: string
            format: input-path
            description: A Visium directory/.h5, a Stereo-seq .gef/.gem, a Visium HD bin table, or an .h5ad the notebook starts from.
          language:
            type: string
            enum: [python]
            description: Notebook language and kernel (python only; the R/OSTA route is a future sibling).
          inject:
            type: array
            items:
              type: string
              format: input-path
            description: Python (.py) scripts to add verbatim as code cells, in order.
          merge:
            type: array
            items:
              type: string
            description: >-
              Notebooks to splice in after the rendered steps, in order. Whole
              notebook (path.ipynb) or a cell range (path.ipynb[0:6]). Must be python.
          config:
            type: string
            format: input-path
            description: JSON file of parameter overrides; keys mirror the parameter names in assets/steps.json.
          output_object:
            type: string
            description: Where the notebook saves its final object, relative to the notebook folder (default spatial.h5ad).
          sample_column:
            type: string
            description: obs column naming each unit's sample (default sample).
          platform:
            type: string
            enum: [auto, visium, visium-hd, stereoseq]
            description: Force a platform (default auto).
          bin_size:
            type: integer
            description: Stereo-seq / Visium HD bin side; null reads the file's bin.
          cellbin:
            type: boolean
            description: Read Stereo-seq cellbin instead of square bins.
          min_genes:
            type: integer
            description: Initial minimum genes per spatial unit (default 100).
          mt_max:
            type: number
            description: Initial maximum mitochondrial percent (default 20).
          mt_pattern:
            type: string
            description: Mitochondrial gene regex (default ^MT- for human; use ^mt- for mouse).
          n_hvg:
            type: integer
            description: Number of variable genes (default 2000).
          n_pcs:
            type: integer
            description: Principal components to compute (default 30).
          dims:
            type: integer
            description: PCs used for neighbours and UMAP (default 30).
          resolution:
            type: number
            description: Clustering resolution (default 0.8).
          seed:
            type: integer
            description: Random seed for reproducible clustering and embedding (default 42).
          neighbor_method:
            type: string
            enum: [knn, radius, delaunay, grid]
            description: Spatial adjacency geometry (default knn).
          n_neighs:
            type: integer
            description: Neighbours for knn/grid adjacency (default 6).
          coord_key:
            type: string
            description: obsm key holding spatial coordinates (default spatial).
          group_by:
            type: string
            description: obs column defining the domains (default leiden).
          svg_method:
            type: string
            enum: [moran, geary]
            description: Spatially variable gene statistic (default moran).
          n_perms:
            type: integer
            description: Permutations for neighbourhood enrichment and ligrec (default 1000).
          reference_path:
            type: string
            format: input-path
            description: Single-cell reference .h5ad for deconvolution; null skips it.
          interactions_path:
            type: string
            format: input-path
            description: Ligand-receptor table for communication; null uses the shipped offline set.
          ccc_threshold:
            type: number
            description: p-value threshold for significant interactions (default 0.01).
          group_column:
            type: string
            description: Condition/sample column for cross-sample differential; null skips it.
          kernel:
            type: string
            description: Kernel name written to the notebook (default plugin-spatial-transcriptomics-spatial).
          self_check:
            type: boolean
            description: Execute every cell once in the output directory and fail if one breaks.
      approval: write
      timeoutSeconds: 1800
    - name: py_env_probe
      description: >-
        Report which spatial Python packages the current environment provides,
        list the packages each workflow module needs, and say which modules can
        run here and which optional sections cannot. Call this before promising
        a workflow such as model-based deconvolution or liana communication. The
        plugin environment is read-only and prebuilt, so a missing package needs
        a project environment, not an in-environment install.
      run: [python, ./scripts/env_probe.py]
      args:
        type: object
        additionalProperties: false
        properties:
          module:
            type: integer
            description: Restrict the report to this module number (1-6).
          write_panel:
            type: string
            format: project-path
            description: Write a starter marker-panel JSON here for the annotate step.
      approval: read
      output: ./schemas/env-probe.json
      timeoutSeconds: 300
---

# spatial transcriptomics workflow

Run a reproducible spatial analysis chain on one AnnData object at a time.
Every step reads an `.h5ad` (or a raw spatial matrix) and writes one, so steps
chain, can be re-run in isolation, and record the parameters they used. This is
the Python route of the spatial plugin (scanpy 1.12 + squidpy 1.8 +
spatialdata); the R/Bioconductor OSTA route is a future sibling skill and does
not exist here. The analysis structure and teaching text are an original
distillation of OSTA (Crowell et al., CC BY 4.0); the code is original to this
plugin and runs in the private `plugin:spatial` environment.

## Before the first step

Call `st_py_inspect` on an unfamiliar input: it reports the platform (Visium /
Visium HD / Stereo-seq), the bin or cellbin, whether spatial coordinates are
present, and which steps the object is ready for. Call `st_py_env_probe` before
promising model-based deconvolution (cell2location, tangram-sc, stereopy) or
liana communication. The plugin environment ships `scanpy`, `anndata`,
`squidpy`, `spatialdata`, `spatialdata-io`, `omnipath`, `leidenalg`,
`python-igraph`, `scikit-misc`, `scikit-learn`, `numba`, `seaborn`, `h5py` and
`ipykernel`. It does **not** ship the pip-only model tools (`cell2location`,
`tangram-sc`, `stereopy`, `liana`, `scvi-tools`). Report a gap as a fact and name
the next step (a project environment via `env_request`); never install packages
inside the read-only environment.

## The chain

1. `st_py_inspect` — what platform/bin is this? Which steps is it ready for?
2. `st_py_qc` — mitochondrial fraction, thresholds, distributions, tissue QC maps.
3. `st_py_preprocess` — normalize, variable genes, PCA, Leiden clusters, UMAP, spatial adjacency graph, cluster maps.
4. `st_py_markers` — marker genes per spatial domain, as tables.
5. `st_py_annotate` — labels from a marker panel or an explicit mapping.
6. `st_py_spatial` — spatially variable genes, neighbourhood enrichment / co-occurrence, cross-sample differential spatial patterns.
7. `st_py_deconvolve` — cell-type proportions against a single-cell reference (scipy NNLS).
8. `st_py_communication` — ligand-receptor communication between domains (squidpy ligrec).

Skip a step only when the object already contains its result; `st_py_inspect`
says which results are present. Steps 6-8 need the clusters and lognorm snapshot
from `st_py_preprocess`; deconvolution also needs a reference, communication a
ligand-receptor resource (shipped offline by default).

## Rendering a notebook

`st_py_notebook` renders the fragment library (`assets/fragments/`) into one
finished notebook in the project directory. It is the same deterministic
renderer as the single-cell skill, driven by this skill's `steps.json`: steps,
parameter values, and extra sources arrive as arguments, and placeholders are
filled from them and from the data path. It does not choose steps or thresholds
and does not interpret the user's request — those decisions are yours (see
Rules).

- `module` renders a preset (1-6); `steps` renders the listed steps plus their
  declared dependencies (`assets/steps.json`).
- Every tunable value lands in the first code cell (the config cell) with a
  source tag; the analysis cells reference those variables.
- `inject` adds a script verbatim as a code cell; `merge` splices in another
  notebook whole or in part (`path.ipynb[0:6]`). Sources must be Python — one
  kernel per notebook.
- `self_check: true` executes every cell in the output directory before
  delivery, with the same relative paths the user will run with.

## Rules

- **The platform and the analysis unit come first.** A Visium spot mixes several
  cells; a Stereo-seq square bin aggregates N×N DNBs; a cellbin is a segmented
  cell. Read `st_py_inspect`'s platform/bin report and say which unit you
  analysed and why. **Bin size is a biological decision, not a default** — large
  bins for tissue domains, cellbin or small bins for cellular questions.
  Segmentation for cellbin is done upstream by SAW; this plugin reads its output.
- **The renderer fills; the agent decides.** Which steps to render, which
  thresholds to pass, which bin, how to edit the marker panel, which adjacency
  geometry — these are judgements you make from the user's question and the data,
  then pass as arguments. A default that reached the notebook untouched is a
  decision you did not make.
- **One reading, both paths.** Pass the same thresholds read from `st_py_qc` to
  the script chain (as flags) and to `st_py_notebook` (as flags or a `config`
  JSON), so the two paths agree by construction, not by luck.
- **Thresholds are arguments, not guesses.** Read `distribution` from `st_py_qc`
  and pass a value the data supports; say which number you chose and why. A
  threshold that removes most units is a warning (`kept_fraction` below half),
  not a result. In spatial data the low shoulder is often off-tissue edge, not
  empty droplets — read the tissue QC map before filtering anatomy away.
- **Counts come first.** QC, HVG (seurat_v3) and deconvolution need the `counts`
  layer. If `counts_present` is false, re-import the raw matrix rather than
  normalizing normalized values twice.
- **The lognorm snapshot is the expression of record.** Marker tests, panel
  scoring, SVG and ligrec read `layers["lognorm"]` (pre-scale, log-normalized);
  never score scaled values.
- **Coordinates and the adjacency graph define every spatial statistic.** SVG,
  neighbourhood enrichment, co-occurrence and communication read
  `obsm[coord_key]` and the spatial graph. Say which `neighbor_method` and
  parameters you used (knn/radius/delaunay/grid, n_neighs, radius); a different
  adjacency definition can change the conclusion. Without coordinates an object
  supports only non-spatial steps.
- **A spatial label is a region, not a pure cell.** After annotation, say how
  many domains got each label and how many scored below `min_score`; an
  unassigned domain is information. Read the dot plot before accepting the call.
- **Deconvolution needs a reference, and is an estimate.** The reference comes
  from the single-cell plugin (shared counts+obs or pseudobulk handoff). State
  its cell-type granularity and source. scipy NNLS is a linear approximation, not
  a validated model; for cell2location use a project environment. A spot's
  composition is estimated, not observed.
- **SVG and communication are screens, not destiny.** Moran's I ranks spatial
  structure; ligrec ranks candidate ligand-receptor pairs by a permutation test.
  Report the statistic, the correction and the adjacency they came from; do not
  claim a mechanism without orthogonal evidence.
- **Cross-sample comparisons need replicates.** The sample is the experimental
  unit, not the spot. Fewer than two biological replicates per condition means no
  condition-level claim; `st_py_spatial`'s differential section is a descriptive
  screen of changing spatial structure, not a significance test of the difference.

## Reading the results

The JSON each tool returns is the audit record: shapes before and after, the
filters that fired with their unit counts, cluster sizes, variance explained,
the spatial graph density, the SVG ranking, enriched domain pairs, deconvolution
proportions, communication pairs, and every note about what could not run. Quote
those numbers in your answer instead of re-deriving them. A section carrying
`{"skipped": reason}` did not run — never read its absence as a clean result.

See [implementation notes](references/maintainers/implementation-notes.md) for
what each step computes and how OSTA maps onto it, and
[osta-map.md](references/maintainers/osta-map.md) for the chapter-by-chapter
correspondence. For the per-stage decision procedure — which tool to call, which
JSON fields and figures to read, how to set thresholds and adjacency, and the red
lines — read [analysis-playbook.md](references/analysis-playbook.md).
[references/README.md](references/README.md) indexes which document serves which
audience; [references/step-notes/](references/step-notes/) carries the per-step
learner-facing short version (the notebook's teaching cells).