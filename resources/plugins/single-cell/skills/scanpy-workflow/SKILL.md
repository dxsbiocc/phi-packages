---
name: scanpy-workflow
description: >-
  Run a scanpy single-cell RNA-seq workflow on an AnnData object, .h5ad, or 10x
  matrix: inspect the input, quality control with scrublet, normalization,
  variable genes, PCA, batch integration (harmony/combat), Leiden clustering,
  UMAP, marker genes, cell-type annotation from a marker panel, pseudobulk
  aggregation, PAGA/diffusion-pseudotime trajectory inference, heteromeric
  ligand-receptor cell-cell communication screening, and AUCell/ULM TF regulon
  activity inference. Use for 单细胞分析, scRNA-seq, scanpy, Python, 质控,
  降维聚类, 批次整合, marker 基因, 细胞类型注释, pseudobulk, 轨迹推断, 细胞通讯,
  配受体, 转录因子活性, regulon, AUCell, ULM, or when a .h5ad / 10x matrix needs
  a reproducible analysis chain. The Python counterpart of seurat-workflow; not
  for spatial transcriptomics or scATAC-seq.
phi:
  environment: plugin:scanpy
  attachTo:
    - SingleCell
  scripts:
    - name: py_inspect
      description: >-
        Read a .h5ad, a 10x Cell Ranger matrix directory or .h5, a .loom, or a
        .csv/.tsv count table and report its shape, layers, embeddings, and the
        metadata columns worth using. Lists which follow-up steps the object
        is ready for. Call this first on an unfamiliar input. With output,
        saves the loaded object as .h5ad so later steps can chain.
      run: [python, ./scripts/project_inspect.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A .h5ad, a 10x matrix directory or .h5, a .loom, or a .csv/.tsv count table.
          output:
            type: string
            format: project-path
            description: Where to save the loaded AnnData object as .h5ad, inside the project.
          min_cells:
            type: integer
            minimum: 0
            description: Drop genes detected in fewer cells.
          min_features:
            type: integer
            minimum: 0
            description: Drop cells with fewer genes.
          var_names:
            type: string
            enum: [gene_symbols, gene_ids]
            description: Which 10x feature column becomes var_names (default gene_symbols).
      approval: read
      output: ./schemas/step-result.json
      timeoutSeconds: 600
    - name: py_qc
      description: >-
        Compute mitochondrial, ribosomal, and hemoglobin fractions, detect
        doublets with scrublet when asked, apply the given thresholds, and
        report how many cells each threshold removed plus the metric
        distributions behind a sensible choice. Writes before/after violin
        plots. Start from raw counts.
      run: [python, ./scripts/qc_filter.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: AnnData saved as .h5ad, with a counts layer.
          output:
            type: string
            format: project-path
            description: Where to save the filtered object as .h5ad.
          figures_dir:
            type: string
            format: project-path
            description: Directory for QC plots; omit to skip figures.
          feature_low:
            type: number
            description: Minimum genes per cell (n_genes_by_counts).
          feature_high:
            type: number
            description: Maximum genes per cell (n_genes_by_counts).
          count_low:
            type: number
            description: Minimum UMI counts per cell (total_counts).
          count_high:
            type: number
            description: Maximum UMI counts per cell (total_counts).
          mt_max:
            type: number
            description: Maximum mitochondrial percent per cell (default 10).
          mt_pattern:
            type: string
            description: Regex for mitochondrial genes (default ^MT- for human; use ^mt- for mouse).
          detect_doublets:
            type: boolean
            description: Run scrublet and drop cells called doublets.
          doublet_expected_rate:
            type: number
            description: Expected doublet rate; default scales with cell number.
          quantiles:
            type: array
            items:
              type: number
            description: Quantiles to report per metric (default 0.01 and 0.99).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 1800
    - name: py_preprocess
      description: >-
        Normalize (log1p of library-size-normalized counts), select variable
        genes on the counts layer, snapshot expression to the lognorm layer,
        scale, run PCA
        on the variable genes, build the neighbor graph, cluster with Leiden at
        one or more resolutions, and embed with UMAP and/or t-SNE. Reports
        variance explained, cluster sizes, and the top variable genes. Writes
        a UMAP and an elbow plot.
      run: [python, ./scripts/preprocess_cluster.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: AnnData saved as .h5ad.
          output:
            type: string
            format: project-path
            description: Where to save the processed object as .h5ad.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the UMAP and elbow plots.
          normalization:
            type: string
            enum: [lognormalize]
            description: Normalization method; SCTransform lives in the Seurat route.
          scale_factor:
            type: number
            description: Library-size target sum (default 10000).
          n_variable_features:
            type: integer
            description: Number of variable genes (default 2000).
          hvg_flavor:
            type: string
            description: HVG flavor on the counts layer (default seurat_v3, the vst analogue).
          n_pcs:
            type: integer
            description: Principal components to compute (default 30).
          dims_use:
            type: integer
            description: Dimensions used for neighbors and embedding (default n_pcs).
          resolution:
            type: array
            items:
              type: number
            description: One or more Leiden resolutions, for example 0.3 0.8.
          regress_vars:
            type: array
            items:
              type: string
            description: Metadata columns to regress out before scaling.
          batch_key:
            type: string
            description: obs column to select variable genes per batch.
          umap:
            type: boolean
            description: Run UMAP (default true).
          tsne:
            type: boolean
            description: Also run t-SNE.
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: py_integrate
      description: >-
        Correct a batch effect with harmony or combat (bbknn is a documented
        gap: no python 3.13 build), store the corrected embedding as
        X_pca_integrated, and recompute neighbors, Leiden clusters and UMAP on
        it. Writes before/after UMAPs colored by batch. Run after
        sc_py_preprocess when a batch split is suspected; downstream tools then
        read corrected clusters.
      run: [python, ./scripts/integrate.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: AnnData saved as .h5ad, after sc_py_preprocess.
          output:
            type: string
            format: project-path
            description: Where to save the integrated object as .h5ad.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the before/after batch UMAPs.
          method:
            type: string
            enum: [harmony, combat, bbknn]
            description: Integration method (default harmony); bbknn reports the environment gap.
          batch_key:
            type: string
            description: obs column naming the batch (default sample).
          resolution:
            type: array
            items:
              type: number
            description: Leiden resolutions for the re-clustering (default 0.8).
          dims_use:
            type: integer
            description: Dimensions used for neighbors (default 30).
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: py_markers
      description: >-
        Rank marker genes per group with rank_genes_groups on the pre-scale
        lognorm snapshot, write one CSV per group, and draw a marker heatmap.
        Reports the top genes and log2 fold changes per group. The same
        threshold flags as the Seurat route filter the tables.
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
            description: Directory for one markers-<group>.csv per cluster.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the marker heatmap.
          group_by:
            type: string
            description: obs column to group by (default leiden).
          test:
            type: string
            enum: [wilcoxon, t-test, t-test_overestim_var, logreg]
            description: Ranking test (default wilcoxon).
          logfc_threshold:
            type: number
            description: Minimum log2 fold change (default 0.25).
          min_pct:
            type: number
            description: Minimum fraction of cells expressing the gene (default 0.1).
          only_positive:
            type: boolean
            description: Report positive markers only (default true).
          top_n:
            type: integer
            description: Top genes to report per group (default 10).
          heatmap:
            type: boolean
            description: Draw the marker heatmap (default true).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: py_annotate
      description: >-
        Score every cluster against a JSON marker panel (cell type -> gene
        symbols) with the same panel-mean-minus-rest rule as the Seurat route,
        and label clusters either from the highest score or from an explicit
        mapping JSON (cluster -> label). Writes the label column into obs and
        draws a marker dot plot. Reports per-cluster scores and which panel
        genes were absent. Optionally runs CellTypist when it is installed.
      run: [python, ./scripts/annotate_clusters.py]
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
            description: Where to save the annotated object as .h5ad.
          panel:
            type: string
            format: input-path
            description: JSON object mapping cell type to an array of gene symbols.
          mapping:
            type: string
            format: input-path
            description: JSON object mapping every cluster to a label; overrides panel scoring.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the marker dot plot.
          group_by:
            type: string
            description: Cluster column (default leiden).
          label_column:
            type: string
            description: obs column to write labels into (default cell_type).
          method:
            type: string
            enum: [panel, celltypist]
            description: Annotation method; celltypist needs the package and a model installed.
          min_score:
            type: number
            description: Minimum panel score to accept a label (default 0).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: py_pseudobulk
      description: >-
        Sum the counts layer per sample (optionally per sample x group) into a
        pseudobulk counts table plus a sample sheet, in the exact column layout
        the repository's differential-expression wrappers read (gene_id first,
        sample-named columns, a sheet whose sample column matches them). The
        result carries wrapper_hint: pass those values to the DESeq2/limma/edgeR
        wrapper to test the contrast. Reports the design, dropped genes, and
        whether enough biological replicates exist.
      run: [python, ./scripts/pseudobulk.py]
      args:
        type: object
        additionalProperties: false
        required: [input, sample_column]
        properties:
          input:
            type: string
            format: input-path
            description: AnnData saved as .h5ad, with a counts layer.
          sample_column:
            type: string
            description: obs column naming the sample of each cell.
          group_column:
            type: string
            description: obs column naming the condition, for a sample-level comparison.
          counts_output:
            type: string
            format: project-path
            description: Path for the tab-separated pseudobulk counts table.
          samples_output:
            type: string
            format: project-path
            description: Path for the tab-separated sample sheet.
          min_cells:
            type: integer
            description: Minimum cells per sample to keep it (default 10).
          min_count:
            type: number
            description: Count a gene as expressed at or above this value (default 1).
          min_cells_expressing:
            type: integer
            description: Keep genes expressed in at least this many units (default 3).
          min_samples_per_group:
            type: integer
            description: Replicates expected per group for a valid comparison (default 2).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 1800
    - name: py_trajectory
      description: >-
        Trajectory inference on a clustered object: PAGA graph of cluster
        relationships plus diffusion pseudotime rooted at a chosen cluster
        (default the first). Needs a neighbor graph, so run after
        sc_py_preprocess or sc_py_integrate. The root is a biological judgement:
        pick the cluster that starts the process, not the first by number.
      run: [python, ./scripts/trajectory.py]
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
            description: Where to save the object with pseudotime attached.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the PAGA and pseudotime UMAP plots.
          group_by:
            type: string
            description: Cluster column (default leiden).
          root_cluster:
            type: string
            description: Cluster the pseudotime starts from (default the first level).
          n_dcs:
            type: integer
            description: Diffusion components for pseudotime (default 10).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 1800
    - name: py_communication
      description: >-
        Infer putative heteromeric ligand-receptor cell-cell communication axes
        between annotated clusters or cell types on the unscaled lognorm layer
        (CellPhoneDB subunit-minimum expression rule + cell-label permutation
        p-values and BH FDR), writing an interactions TSV, a sender-receiver hub
        heatmap, and a top L-R dotplot. Uses a curated starter catalog unless an
        OmniPath/CellPhoneDB TSV/CSV is supplied via interactions_path.
      run: [python, ./scripts/communication.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: AnnData saved as .h5ad, after clustering/annotation.
          output:
            type: string
            format: project-path
            description: Optional .h5ad path to save the object with communication metadata attached.
          interactions_output:
            type: string
            format: project-path
            description: Optional TSV path for all active tested ligand-receptor interactions.
          figures_dir:
            type: string
            format: project-path
            description: Directory for communication-heatmap.png and communication-dotplot.png.
          group_by:
            type: string
            description: obs column defining sender/receiver populations (default cell_type).
          interactions_path:
            type: string
            format: input-path
            description: Optional custom ligand-receptor TSV/CSV (ligand, receptor or source, target; '+' for heteromers).
          min_pct:
            type: number
            description: Minimum expressed cell fraction per subunit in sender and receiver groups (default 0.10).
          n_perms:
            type: integer
            minimum: 50
            description: Cell-label permutations for empirical p-value estimation (default 500).
          fdr:
            type: number
            description: BH FDR cutoff for significant interactions (default 0.05).
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 1800
    - name: py_regulon
      description: >-
        Infer per-cell transcription factor (TF) regulon activity using AUCell
        (normalized recovery AUC over top expressed genes) or weighted ULM linear
        modeling on the unscaled lognorm layer, storing X_regulon in obsm and
        computing per-group Mann-Whitney U differential activity (BH FDR) plus a
        top-TF group heatmap.
      run: [python, ./scripts/regulon.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: Clustered or annotated AnnData object saved as .h5ad.
          output:
            type: string
            format: project-path
            description: Where to save the updated AnnData object with obsm["X_regulon"].
          regulon_output:
            type: string
            format: project-path
            description: Output TSV path for cell-by-TF activity scores (default regulon/regulon_activity.tsv).
          figures_dir:
            type: string
            format: project-path
            description: Directory for TF regulon heatmap PNG; omit to skip figure.
          group_by:
            type: string
            description: obs column grouping cells for differential regulon summary (default leiden).
          regulons_path:
            type: string
            format: input-path
            description: Optional TSV/CSV (source/tf, target/gene, weight) or .gmt; omit for curated 18-TF starter panel.
          regulon_method:
            type: string
            enum: [aucell, ulm]
            description: Regulon activity scoring algorithm (default aucell).
          top_frac:
            type: number
            description: Top expressed gene fraction for AUCell recovery curve (default 0.05).
          min_targets:
            type: integer
            minimum: 1
            description: Minimum matched target genes required per TF (default 5).
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 1800
    - name: py_notebook
      description: >-
        Render a runnable analysis notebook into the project from the fragment
        library, then open it with Phi's notebook tools and run cells on demand.
        The template is rendered, not copied: placeholders are filled from the
        input data and the requested steps, and every tunable value lands in one
        config cell at the top that the analysis cells reference. Pick a preset
        with module (1 QC, 2 plus clustering, 3 plus integration, 4 plus
        markers, 5 plus pseudobulk, 6 plus annotation, 7 plus trajectory, 8 plus
        communication, 9 plus regulon) or
        name steps directly; dependencies are added automatically. inject adds
        an external Python or R script as a code cell; merge splices in another
        notebook, whole or a cell range. All sources must match the notebook
        language (one kernel per notebook). Starts from a count matrix, .h5, or
        .h5ad: raw-read processing (alignment, Cell Ranger) is the wrapper
        tree's job, not a notebook cell. Use self_check to execute every cell
        in the output directory now and fail loudly if one breaks.
      run: [python, ./scripts/build_notebook.py]
      args:
        type: object
        additionalProperties: false
        required: [output, data]
        properties:
          module:
            type: integer
            enum: [1, 2, 3, 4, 5, 6, 7, 8, 9]
            description: NBIS/extended lab preset (1-9). Pass module or steps, not both.
          steps:
            type: array
            items:
              type: string
              enum: [setup, load, qc-metrics, qc-plots, doublets, filter, normalize, hvg, scale, pca, integrate, batch-compare, cluster, markers, annotate, assign-labels, trajectory, communication, regulon, pseudobulk, save]
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
            description: A 10x matrix directory, .h5, or .h5ad the notebook starts from.
          language:
            type: string
            enum: [python, r]
            description: >-
              Notebook language and kernel (default python). The Python fragment
              library backs module/steps; an r notebook is assembled from inject
              and merge sources only.
          inject:
            type: array
            items:
              type: string
              format: input-path
            description: Python (.py) or R (.R) scripts to add verbatim as code cells, in order.
          merge:
            type: array
            items:
              type: string
            description: >-
              Notebooks to splice in after the rendered steps, in order. Whole
              notebook (path.ipynb) or a cell range (path.ipynb[0:6], indices
              over all cells). Must match the notebook language.
          config:
            type: string
            format: input-path
            description: JSON file of parameter overrides; keys mirror the parameter names below.
          output_object:
            type: string
            description: Where the notebook saves its final object, relative to the notebook folder (default analyzed.h5ad).
          sample_column:
            type: string
            description: obs column naming each cell's sample (default sample).
          min_genes:
            type: integer
            description: Initial minimum genes per cell (default 200; lower it for shallow data).
          max_counts:
            type: integer
            description: Initial maximum UMI counts per cell (default 25000).
          mt_max:
            type: number
            description: Initial maximum mitochondrial percent (default 10).
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
            description: Dimensions used for neighbours and UMAP (default 30).
          resolution:
            type: number
            description: Clustering resolution (default 0.8).
          seed:
            type: integer
            description: Random seed for reproducible clustering and embedding (default 42).
          batch_key:
            type: string
            description: obs column naming the batch; integration skips when absent (default sample).
          batch_method:
            type: string
            enum: [harmony, combat]
            description: Integration method in the notebook (default harmony).
          hvg_flavor:
            type: string
            description: HVG flavor on the counts layer (default seurat_v3).
          kernel:
            type: string
            description: Kernel name written to the notebook (default plugin-single-cell-scanpy).
          self_check:
            type: boolean
            description: Execute every cell once in the output directory and fail if one breaks.
      approval: write
      output: ./schemas/notebook-result.json
      timeoutSeconds: 1800
    - name: py_env_probe
      description: >-
        Report which single-cell Python packages the current environment
        provides, list the packages each NBIS lab needs, and say which lab
        cores can run here and which optional sections cannot. Call this before
        promising a workflow such as scVI integration or spatial analysis. The
        plugin environment is read-only and prebuilt, so a missing package
        needs a project environment, not an in-environment install.
      run: [python, ./scripts/env_probe.py]
      args:
        type: object
        additionalProperties: false
        properties:
          module:
            type: integer
            description: Restrict the report to this lab number (1-7).
          write_panel:
            type: string
            format: project-path
            description: Write a starter marker-panel JSON here for the annotate step.
      approval: read
      output: ./schemas/env-probe.json
      timeoutSeconds: 300
---

# scanpy single-cell workflow

Run a reproducible scanpy analysis chain on one AnnData object at a time.
Every step reads an `.h5ad` and writes one, so steps chain, can be re-run in
isolation, and record the parameters they used. This is the Python route of
the single-cell plugin; the sibling `seurat-workflow` skill runs the same
analysis on Seurat objects. Route by the user's language preference and input
format (see the SingleCell agent rules): Python / `.h5ad` here, R / `.rds`
there. The two routes share threshold semantics, panel scoring, and the
pseudobulk handoff layout, so a project can switch routes without rewriting
its decisions.

## Before the first step

Call `sc_py_env_probe` on an unfamiliar machine, or whenever a requested
workflow is scVI integration, CellTypist annotation, enrichment, or spatial
analysis. It reports what the environment actually provides. The plugin
environment (`plugin:scanpy`, Python 3.13 + scanpy 1.13) ships `scanpy`,
`anndata`, `leidenalg`, `python-igraph`, `harmonypy`, `scrublet`,
`scikit-misc`, `scikit-learn`, `seaborn`, and `ipykernel`. It does **not** ship
`bbknn` (no python 3.13 build on bioconda), `celltypist`, `scanorama`,
`scvi-tools`, `gseapy`, or `squidpy`. Report a gap as a fact and name the next
step; never install packages inside the read-only environment.

## The chain

1. `sc_py_inspect` — what is in the file? Which steps is it ready for?
2. `sc_py_qc` — scrublet doublets, mitochondrial fraction, thresholds, distributions.
3. `sc_py_preprocess` — normalize, variable genes, PCA, Leiden clusters, UMAP.
4. `sc_py_integrate` — harmony/combat batch correction, then re-cluster (only when a batch split is suspected).
5. `sc_py_markers` — marker genes per cluster, as tables.
6. `sc_py_annotate` — labels from a marker panel or an explicit mapping.
7. `sc_py_pseudobulk` — sample-level counts for differential expression.
8. `sc_py_trajectory` — PAGA and diffusion pseudotime (when the question is a process, not groups).
9. `sc_py_communication` — heteromeric ligand-receptor permutation screening across cell types.
10. `sc_py_regulon` — AUCell or ULM TF regulon activity inference across single cells.

Skip a step only when the object already contains its result; `sc_py_inspect`
says which results are present.

## Rendering a notebook

`sc_py_notebook` renders the fragment library (`assets/fragments/`) into one
finished notebook in the project directory. It is the same deterministic
renderer as the Seurat skill (byte-identical `build_notebook.py`), driven by
this skill's `steps.json`: steps, parameter values, and extra sources arrive as
arguments, and placeholders are filled from them and from the data path. It
does not choose steps or thresholds and does not interpret the user's request
— those decisions are yours (see Rules).

- `module` renders an NBIS lab preset (1-7); `steps` renders the listed steps
  plus their declared dependencies (`assets/steps.json`).
- Every tunable value lands in the first code cell (the config cell) with a
  source tag; the analysis cells reference those variables.
- `inject` adds a script verbatim as a code cell; `merge` splices in another
  notebook whole or in part (`path.ipynb[0:6]`). All sources must match
  `language` — one kernel per notebook.
- `self_check: true` executes every cell in the output directory before
  delivery, with the same relative paths the user will run with.

## Rules

- **R or Python is the user's route, not yours.** This skill is the Python/scanpy
  half of the plugin; `seurat-workflow` (`sc_*`, `plugin:seurat`) is the R half.
  Follow the user's stated preference; default to Python for `.h5ad` input and
  ask when only a 10x matrix is given. Never mix routes in one analysis — hand
  off only through pseudobulk tables or the 10x-style exchange in
  `references/cross-route-handoff.md`.
- **The renderer fills; the agent decides.** Which steps to render, which
  thresholds to pass, how to edit the marker panel, which cluster roots a
  trajectory — these are judgements you make from the user's question and the
  data, then pass as arguments. The renderer never infers a threshold or picks
  a step from the data; a default that reached the notebook untouched is a
  decision you did not make.
- **One reading, both paths.** Pass the same thresholds read from `sc_py_qc`
  to the script chain (as flags) and to `sc_py_notebook` (as flags or a
  `config` JSON), so the two paths agree by construction, not by luck. The
  same numbers apply when the project switches to the Seurat route.
- **Thresholds are arguments, not guesses.** The labs choose cutoffs by eye.
  Read `distribution` from `sc_py_qc` and pass a value the data supports; say
  which number you chose and why.
- **A threshold that removes most cells is a warning, not a result.** The QC
  result reports `kept_fraction`; below half, stop and reconsider.
- **Counts come first.** QC, HVG (seurat_v3), and pseudobulk need the `counts`
  layer. If `counts_present` is false, re-import the raw matrix rather than
  normalizing normalized values twice.
- **The lognorm snapshot is the expression of record.** Marker tests and panel
  scoring read `layers["lognorm"]` (pre-scale, log-normalized); never score
  scaled values.
- **Integrate before interpreting batches.** With more than one sample, check
  the batch UMAP from `sc_py_preprocess`; when batches separate, run
  `sc_py_integrate` and read markers from the corrected clusters. Judge
  integration by the before/after plots: mixing without erasing biology.
- **Report the design, not just the label.** After annotation, say how many
  clusters got each label and how many scored below `min_score`; an unassigned
  cluster is information.
- **Pseudotime is relative, not destiny.** Name the root cluster and why;
  report unreachable cells instead of hiding them.
- **`avg_log2FC` is not a biological effect size.** Pass marker tables to
  downstream interpretation; do not claim a fold change without the test and
  the comparison it came from.

## Reading the results

The JSON each tool returns is the audit record: shapes before and after, the
filters that fired with their cell counts, cluster sizes, variance explained,
and every note about what could not run. Quote those numbers in your answer
instead of re-deriving them.

See [implementation notes](references/maintainers/implementation-notes.md) for
what each step computes and how the NBIS labs map onto it, and
[nbis-porting.md](references/maintainers/nbis-porting.md) for the decisions
behind this port. [cross-route-handoff.md](references/cross-route-handoff.md)
covers moving objects between the two routes. For the per-stage decision
procedure — which tool to call, which JSON fields and figures to read, how to
set thresholds, and the red lines — read
[analysis-playbook.md](references/analysis-playbook.md).
[references/README.md](references/README.md) indexes which document serves
which audience; [references/step-notes/](references/step-notes/) carries the
per-step learner-facing short version (the notebook's teaching cells).