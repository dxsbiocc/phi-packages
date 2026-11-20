# Pipeline notes

What each `sc_*` tool computes, and which scNotebooks cell it came from.
Read this when adapting a step rather than re-running it. Since 0.6.0 the
analysis lives in `assets/fragments/` (per-step explanations in
`references/step-notes/`) — the tools are thin adapters over the same
fragments the notebook renders, so adapt the fragment, not the script.

## `sc_inspect` — what is in this file?

Source: Module 04 "Data input, sparse UMI-count matrix", Module 09 "Read in raw
gene counts and metadata".

The notebooks open a file with `Read10X_h5()` and then interrogate it by hand
(`dim`, `class`, `rownames`, `colnames`). The tool does that interrogation once
and reports:

- shape, assays, default assay, layers, reductions, graphs;
- whether a `counts` layer exists — without it QC and pseudobulk are invalid;
- metadata columns likely to matter downstream (`orig.ident`, `nCount_RNA`,
  `nFeature_RNA`, `percent.mt`, cluster and annotation columns, `sample`,
  `condition`, `batch`), with numeric ranges or category counts;
- `candidate_steps`, the steps this object is already capable of.

`.h5` and `.h5ad` are refused with an explanation rather than a crash: reading
`.h5` needs `hdf5r`, which the plugin environment does not ship, and `.h5ad` is AnnData, which
belongs to the `scanpy` skill.

## `sc_qc` — what does a defensible threshold look like?

Source: Module 04 "Quality control (QC)", "Identifying Doublets", "Ambient RNA
Correction (Optional)".

The notebook:

1. filters genes and cells at object creation (`min.cells = 3`,
   `min.features = 200`);
2. runs `scDblFinder` on a `SingleCellExperiment` and keeps singlets;
3. adds `percent.mt` with `PercentageFeatureSet(sampleA, pattern = "^MT-")`;
4. plots violins, histograms, and feature scatter plots;
5. filters with `nCount_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 10`.

The tool keeps the metrics and the filters but makes the numbers explicit and
auditable. It reports the distribution (`min`, `q01`, `median`, `q99`, `max`) of
each metric so a threshold can be chosen from the data, then reports, per
filter, how many cells it removed, and the surviving fraction. A run that keeps
under half the cells produces a warning instead of quietly proceeding.

Mitochondrial/ribosomal/hemoglobin patterns are arguments, because human `^MT-`
and mouse `^mt-` differ and the notebooks hard-code the human form.

Doublet detection is best effort: `scDblFinder` when installed (with a rate
derived from cell number unless given), otherwise a note saying it did not run.
The tool never claims doublets were removed when they were not.

## `sc_preprocess` — one pass instead of eight cells

Source: Module 04 "Normalization", "Identify variable genes", "Scaling gene
expression", "PCA", "UMAP", "Clustering", and the SCTransform section.

The notebooks walk these as separate cells so a learner can inspect each. A
pipeline needs them in one deterministic pass, so the tool chains:
`NormalizeData` → `FindVariableFeatures` → `ScaleData` → `RunPCA` →
`FindNeighbors` → `FindClusters` → `RunUMAP`.

`--normalization sct` swaps the first three for `SCTransform`, matching the
notebook's "Alternative Normalization" section.

Reported: variance explained per PC and for the first ten (the elbow heuristic),
cluster sizes, top variable features, and every parameter used. Repeated
`--resolution` values produce `seurat_clusters_<value>` columns, so a resolution
scan is one call instead of a loop.

## `sc_markers` — tables, not just a heatmap

Source: Module 04 "Differential Expression", Module 10 marker dot plots, Module
11 marker peaks (for RNA).

`FindAllMarkers` with the notebook's thresholds, then one CSV per cluster ordered
by adjusted p-value then descending log2 fold change. The full table is stored on
the object (`misc$markers`) so `sc_annotate` and any later tool can read it
without recomputing.

The heatmap is illustrative; the tables are the result. `avg_log2FC` is reported
with the test name and the comparison so it is not mistaken for an effect size.

## `sc_annotate` — panel scoring or an explicit mapping

Source: Module 04 "Annotating Cell Types" (automated: SingleR + celldex; manual:
a named marker list and a dot plot).

The notebooks keep the panel in a code cell. The tool takes a JSON panel
(`{"T_cell": ["CD3D", ...]}`) and scores every cluster as *mean expression of
the panel genes minus the mean of all other genes in that cluster*. Using a
difference of means keeps panels of different sizes comparable, which a plain sum
does not.

Two ways to finish:

- **panel only** — each cluster gets its highest-scoring cell type, or
  `unassigned` when that score is below `--min-score`;
- **`--mapping`** — an explicit `{"0": "T_cell", ...}`. Every cluster must be
  labelled; a partial mapping is an error, because silently leaving clusters
  unlabelled is how mislabelled data reaches a figure.

`--method singler` runs SingleR when `SingleR` and `celldex` are installed;
otherwise it returns the install gap rather than failing the run.

Panel genes absent from the object are counted and reported — usually a species
naming difference (human `MT-` vs mouse `mt-`).

## `sc_pseudobulk` — the step the notebooks never take

The notebooks analyse one sample and describe cell types; they never compare
conditions at the sample level, which is what a differential-expression test
needs. This step is the bridge to the repository's `differential-expression`
wrappers:

- sum counts per sample (or per sample × group, so a sample split across
  conditions is never merged);
- drop units with too few cells, and genes below an expression floor;
- write a counts table and a sample sheet whose columns match what the wrappers
  expect;
- report the design, including whether each group has enough replicates.

A single sample per group is reported as a warning, because a two-group test
without replicates cannot separate condition from sample.
