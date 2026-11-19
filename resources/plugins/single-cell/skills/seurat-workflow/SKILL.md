---
name: seurat-workflow
description: >-
  Run a Seurat single-cell RNA-seq workflow on a Seurat object or 10x matrix:
  inspect the input, quality control and doublet detection, normalization,
  variable features, PCA, clustering, UMAP/t-SNE, marker genes, cell-type
  annotation from a marker panel, and pseudobulk aggregation for sample-level
  comparison. Use for 单细胞分析, scRNA-seq, Seurat, 质控, 降维聚类, marker 基因,
  细胞类型注释, pseudobulk, or when a .rds / 10x matrix needs a reproducible
  analysis chain. Not for spatial transcriptomics, scATAC-seq, trajectory
  inference, or cell-cell communication (those need packages the plugin
  environment omits; run env_probe first and see the module notes). Also indexes
  the scNotebooks teaching curriculum: use for 单细胞课程, scNotebooks, 教学,
  培训材料, 模块 1-13, 学习路径 (references/maintainers/scnotebooks-curriculum.md).
phi:
  environment: plugin:seurat
  attachTo:
    - SingleCell
  scripts:
    - name: inspect
      description: >-
        Read a Seurat .rds, a 10x Cell Ranger matrix directory, or a
        .csv/.tsv count table and report its shape, assays, layers, reductions,
        and the metadata columns worth using. Lists which follow-up steps the
        object is ready for. Call this first on an unfamiliar input. With
        output, saves the loaded object as .rds so later steps can chain.
      run: [Rscript, ./scripts/project_inspect.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A Seurat .rds, a 10x matrix directory, or a .csv/.tsv count table.
          output:
            type: string
            format: project-path
            description: Where to save the loaded Seurat object as .rds, inside the project.
          project:
            type: string
            description: Project name recorded on a newly created Seurat object.
          min_cells:
            type: integer
            minimum: 0
            description: Drop genes detected in fewer cells (CreateSeuratObject min.cells).
          min_features:
            type: integer
            minimum: 0
            description: Drop cells with fewer genes (CreateSeuratObject min.features).
      approval: read
      timeoutSeconds: 600
    - name: qc
      description: >-
        Compute mitochondrial, ribosomal, and hemoglobin fractions, detect
        doublets when scDblFinder is available, apply the given thresholds, and
        report how many cells each threshold removed plus the metric
        distributions behind a sensible choice. Writes before/after violin
        plots. Start from raw counts.
      run: [Rscript, ./scripts/qc_filter.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: Seurat object saved as .rds, with a counts layer.
          output:
            type: string
            format: project-path
            description: Where to save the filtered object as .rds.
          figures_dir:
            type: string
            format: project-path
            description: Directory for QC plots; omit to skip figures.
          feature_low:
            type: number
            description: Minimum genes per cell (nFeature_RNA).
          feature_high:
            type: number
            description: Maximum genes per cell (nFeature_RNA).
          count_low:
            type: number
            description: Minimum UMI counts per cell (nCount_RNA).
          count_high:
            type: number
            description: Maximum UMI counts per cell (nCount_RNA).
          mt_max:
            type: number
            description: Maximum mitochondrial percent per cell (default 10).
          mt_pattern:
            type: string
            description: Regex for mitochondrial genes (default ^MT- for human; use ^mt- for mouse).
          detect_doublets:
            type: boolean
            description: Run scDblFinder and drop cells called doublets.
          doublet_expected_rate:
            type: number
            description: Expected doublet rate; default scales with cell number.
          quantiles:
            type: array
            items:
              type: number
            description: Quantiles to report per metric (default 0.01 and 0.99).
      approval: write
      timeoutSeconds: 1800
    - name: preprocess
      description: >-
        Normalize (LogNormalize or SCTransform), select variable features,
        scale, run PCA, build the neighbor graph, cluster at one or more
        resolutions, and embed with UMAP and/or t-SNE. Reports variance
        explained, cluster sizes, and the top variable genes. Writes a UMAP and
        an elbow plot.
      run: [Rscript, ./scripts/preprocess_cluster.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: Seurat object saved as .rds.
          output:
            type: string
            format: project-path
            description: Where to save the processed object as .rds.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the UMAP and elbow plots.
          normalization:
            type: string
            enum: [lognormalize, sct]
            description: Normalization method (default lognormalize).
          scale_factor:
            type: number
            description: LogNormalize scale factor (default 10000).
          n_variable_features:
            type: integer
            description: Number of variable features (default 2000).
          selection_method:
            type: string
            description: Variable-feature method, for example vst or mean.var.plot.
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
            description: One or more clustering resolutions, for example 0.3 0.8.
          regress_vars:
            type: array
            items:
              type: string
            description: Metadata columns to regress out while scaling.
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
      timeoutSeconds: 3600
    - name: markers
      description: >-
        Rank marker genes per group with FindAllMarkers, write one CSV per
        group, and draw a marker heatmap. The full table is stored on the
        object so the annotate step can reuse it. Reports the top genes and
        log2 fold changes per group.
      run: [Rscript, ./scripts/find_markers.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: Seurat object saved as .rds, after clustering.
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
            description: Metadata column to group by (default seurat_clusters).
          test:
            type: string
            description: Test passed to FindAllMarkers (default wilcox).
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
          assay:
            type: string
            description: Assay name; default the object's active assay.
      approval: write
      timeoutSeconds: 3600
    - name: annotate
      description: >-
        Score every cluster against a JSON marker panel (cell type -> gene
        symbols) and label clusters either from the highest score or from an
        explicit mapping JSON (cluster -> label). Writes the label column into
        the object and draws a marker dot plot. Reports per-cluster scores and
        which panel genes were absent. Optionally runs SingleR when it is
        installed.
      run: [Rscript, ./scripts/annotate_clusters.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: Seurat object saved as .rds, after clustering.
          output:
            type: string
            format: project-path
            description: Where to save the annotated object as .rds.
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
            description: Cluster column (default seurat_clusters).
          label_column:
            type: string
            description: Metadata column to write labels into (default cell_type).
          method:
            type: string
            enum: [panel, singler]
            description: Annotation method; singler needs SingleR and celldex installed.
          min_score:
            type: number
            description: Minimum panel score to accept a label (default 0).
      approval: write
      timeoutSeconds: 3600
    - name: pseudobulk
      description: >-
        Sum counts per sample (optionally per sample x group) into a pseudobulk
        counts table plus a sample sheet, in the exact column layout the
        repository's differential-expression wrappers read (gene_id first,
        sample-named columns, a sheet whose sample column matches them). The
        result carries wrapper_hint: pass those values to the DESeq2/limma/edgeR
        wrapper to test the contrast. Reports the design, dropped genes, and
        whether enough biological replicates exist.
      run: [Rscript, ./scripts/pseudobulk.R]
      args:
        type: object
        additionalProperties: false
        required: [input, sample_column]
        properties:
          input:
            type: string
            format: input-path
            description: Seurat object saved as .rds, with a counts layer.
          sample_column:
            type: string
            description: Metadata column naming the sample of each cell.
          group_column:
            type: string
            description: Metadata column naming the condition, for a sample-level comparison.
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
      timeoutSeconds: 1800
    - name: notebook
      description: >-
        Render a runnable analysis notebook into the project from the fragment
        library, then open it with Phi's notebook tools and run cells on demand.
        The template is rendered, not copied: placeholders are filled from the
        input data and the requested steps, and every tunable value lands in one
        config cell at the top that the analysis cells reference. Pick a preset
        with module (4 QC-to-pseudobulk, 5 plus batch comparison, 8 plus protein
        assay, 9 plus clonotype join) or name steps directly; dependencies are
        added automatically. inject adds an external R or Python script as a
        code cell; merge splices in another notebook, whole or a cell range.
        All sources must match the notebook language (one kernel per notebook).
        Starts from a count matrix or a saved .rds: raw-read processing
        (alignment, Cell Ranger) is the wrapper tree's job, not a notebook cell.
        Use self_check to execute every cell in the output directory now and
        fail loudly if one breaks.
      run: [python, ./scripts/build_notebook.py]
      args:
        type: object
        additionalProperties: false
        required: [output, data]
        properties:
          module:
            type: integer
            enum: [4, 5, 8, 9]
            description: scNotebooks module preset. Pass module or steps, not both.
          steps:
            type: array
            items:
              type: string
              enum: [setup, load, qc-metrics, qc-plots, doublets, filter, normalize, hvg, scale, pca, cluster, markers, annotate, assign-labels, pseudobulk, save, batch-compare, adt-assay, clonotypes]
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
            description: A 10x matrix directory or a Seurat .rds the notebook starts from.
          language:
            type: string
            enum: [r, python]
            description: >-
              Notebook language and kernel (default r). The R fragment library
              backs module/steps; a python notebook is assembled from inject and
              merge sources only.
          inject:
            type: array
            items:
              type: string
              format: input-path
            description: R (.R) or Python (.py) scripts to add verbatim as code cells, in order.
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
          project:
            type: string
            description: Seurat project name; defaults to the data folder name.
          output_object:
            type: string
            description: Where the notebook saves its final object, relative to the notebook folder (default analyzed.rds).
          sample_column:
            type: string
            description: Metadata column naming each cell's sample (default sample).
          clonotype:
            type: string
            description: Module 9 only - the Cell Ranger contig annotations CSV.
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
          kernel:
            type: string
            description: Kernel name written to the notebook (default plugin-single-cell-seurat / phi-python).
          self_check:
            type: boolean
            description: Execute every cell once in the output directory and fail if one breaks.
      approval: write
      timeoutSeconds: 1800
    - name: env_probe
      description: >-
        Report which single-cell R packages the current environment provides,
        list the packages each scNotebooks module needs, and say which module
        cores can run here and which optional sections cannot. Call this before
        promising a workflow such as trajectory or cell-cell communication. The
        plugin environment is read-only and prebuilt, so a missing package needs
        a project environment, not an in-environment install.
      run: [Rscript, ./scripts/env_probe.R]
      args:
        type: object
        additionalProperties: false
        properties:
          module:
            type: array
            items:
              type: integer
            description: Restrict the report to these module numbers (1-13).
          write_panel:
            type: string
            format: project-path
            description: Write a starter marker-panel JSON here for the annotate step.
      approval: read
      output: ./schemas/env-probe.json
      timeoutSeconds: 300
---

# Seurat single-cell workflow

Run a reproducible Seurat analysis chain on one object at a time. Every step
reads an `.rds` and writes one, so steps chain, can be re-run in isolation, and
record the parameters they used.

## Before the first step

Call `sc_env_probe` on an unfamiliar machine, or whenever a requested
workflow is trajectory, cell-cell communication, spatial, or scATAC-seq. It
reports what the environment actually provides. The plugin environment
(`plugin:seurat`, R 4.5 + Seurat 6) ships `Seurat`, `SeuratObject`,
`SingleCellExperiment`, `jsonlite`, `ggplot2`, `patchwork`, and `ComplexHeatmap`
(plus their dependencies, such as `Matrix` and `uwot`). It does **not** ship
`monocle3`, `liana`, `SoupX`, `DropletUtils`, `semla`, `ArchR`, `scDblFinder`,
`harmony`, or `pheatmap`. Report a gap as a fact and name the next step; never
install packages inside the read-only environment.

## The chain

1. `sc_inspect` — what is in the file? Which steps is it ready for?
2. `sc_qc` — doublets, mitochondrial fraction, thresholds, distributions.
3. `sc_preprocess` — normalize, variable features, PCA, clusters, UMAP.
4. `sc_markers` — marker genes per cluster, as tables.
5. `sc_annotate` — labels from a marker panel or an explicit mapping.
6. `sc_pseudobulk` — sample-level counts for differential expression.

Skip a step only when the object already contains its result; `sc_inspect`
says which results are present.

## Rendering a notebook

`sc_notebook` renders the fragment library (`assets/fragments/`) into one
finished notebook in the project directory. It is a deterministic renderer:
steps, parameter values, and extra sources arrive as arguments, and placeholders
are filled from them and from the data path. It does not choose steps or
thresholds and does not interpret the user's request — those decisions are
yours (see Rules).

- `module` renders a preset; `steps` renders the listed steps plus their
  declared dependencies (`assets/steps.json`).
- Every tunable value lands in the first code cell (the config cell) with a
  source tag; the analysis cells reference those variables.
- `inject` adds a script verbatim as a code cell; `merge` splices in another
  notebook whole or in part (`path.ipynb[0:6]`). All sources must match
  `language` — one kernel per notebook.
- `self_check: true` executes every cell in the output directory before
  delivery, with the same relative paths the user will run with.

## Rules

- **R or Python is the user's route, not yours.** This skill is the R/Seurat
  half of the plugin; `scanpy-workflow` (`sc_py_*`, `plugin:scanpy`) is the
  Python half. Follow the user's stated preference; default to R for `.rds`
  input and ask when only a 10x matrix is given. Never mix routes in one
  analysis — hand off only through pseudobulk tables or the 10x-style exchange
  in scanpy-workflow's `references/cross-route-handoff.md`.
- **The renderer fills; the agent decides.** Which steps to render, which
  thresholds to pass, how to edit the marker panel — these are judgements you
  make from the user's question and the data, then pass as arguments. The
  renderer never infers a threshold or picks a step from the data; a default
  that reached the notebook untouched is a decision you did not make.
- **One reading, both paths.** Pass the same thresholds read from `sc_qc`
  to the script chain (as flags) and to `sc_notebook` (as flags or a
  `config` JSON), so the two paths agree by construction, not by luck.
- **Thresholds are arguments, not guesses.** The notebooks hard-code cutoffs
  chosen by eye. Read `distribution` from `sc_qc` and pass a value that the
  data supports; say which number you chose and why.
- **A threshold that removes most cells is a warning, not a result.** The QC
  result reports `kept_fraction`; below half, stop and reconsider.
- **Counts come first.** QC and pseudobulk need a `counts` layer. If
  `counts_present` is false, re-import the raw matrix rather than normalizing
  normalized values twice.
- **Report the design, not just the label.** After annotation, say how many
  clusters got each label and how many scored below `min_score`; an unassigned
  cluster is information.
- **`avg_log2FC` is not a biological effect size.** Pass marker tables to
  downstream interpretation; do not claim a fold change without the test and the
  comparison it came from.

## Reading the results

The JSON each tool returns is the audit record: shapes before and after, the
filters that fired with their cell counts, cluster sizes, variance explained,
and every note about what could not run. Quote those numbers in your answer
instead of re-deriving them.

See [implementation notes](references/maintainers/implementation-notes.md) for
what each step computes and how the scNotebooks module maps onto it, and
[scnotebooks-porting.md](references/maintainers/scnotebooks-porting.md) for the
decisions behind this port. The teaching curriculum itself is indexed in
[scnotebooks-curriculum.md](references/maintainers/scnotebooks-curriculum.md)
with per-module package facts in
[scnotebooks-module-map.md](references/maintainers/scnotebooks-module-map.md).
For the per-stage decision procedure — which tool to call, which JSON fields
and figures to read, how to set thresholds, and the red lines — read
[analysis-playbook.md](references/analysis-playbook.md).
[references/README.md](references/README.md) indexes which document serves
which audience; [references/step-notes/](references/step-notes/) carries the
per-step learner-facing short version (the notebook's teaching cells).
