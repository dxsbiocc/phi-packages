# Pipeline notes — scanpy route

What each tool computes, and which NBIS workshop-scRNAseq lab it backs. The
chain is deterministic: every step reads one `.h5ad`, writes one, and reports
the parameters it used. Since 0.6.0 the analysis lives in
`assets/fragments/` (per-step explanations in `references/step-notes/`) —
the tools are thin adapters over the same fragments the notebook renders,
so adapt the fragment, not the script.

## py_inspect

Dispatches on the input kind (`.h5ad`, 10x directory/`.h5`, `.loom`, dense
table), snapshots `X` as the `counts` layer when absent, and reports shape,
layers, embeddings, `uns` keys, and the metadata columns worth using.
`candidate_steps` is computed from what the object already contains, so a
re-run skips finished work. A bare `.mtx` is refused: without feature/barcode
files the identifiers would be lost.

## py_qc

`calculate_qc_metrics` on the counts layer with mitochondrial/ribosomal/
hemoglobin gene sets; scrublet doublet scores when `detect_doublets`; then the
threshold filters, each reported with its removed-cell count. The
`distribution` block (min/q01/median/q99/max) is the evidence for choosing
thresholds — the labs eyeball the same violins. `kept_fraction` below 0.5 is a
warning to reconsider, not a result to deliver.

## py_preprocess

`normalize_total` + `log1p`; `highly_variable_genes` with `flavor=seurat_v3`
on the **counts** layer (the analogue of Seurat's vst; needs scikit-misc,
which the environment ships); `layers["lognorm"]` snapshot **before** scaling;
`scale(max_value=10)`; `pca(mask="var.highly_variable")` — the object is never
gene-subset, so downstream steps can still see every gene. Leiden at each
requested resolution (`leiden`, or `leiden_<res>` for several), UMAP/t-SNE.
Variance explained and the elbow plot justify `n_pcs`.

## py_integrate

harmony (`sc.pp.harmony_integrate`) adjusts the PCA embedding;
combat runs on the pre-scale lognorm snapshot and recomputes PCA. Both store
`X_pca_integrated`; neighbours/Leiden/UMAP are then recomputed on it, so
markers and annotation read corrected clusters. bbknn is a documented gap
(no python 3.13 build on bioconda). Before/after UMAPs colored by batch are
the judgement aid: integration should mix samples without erasing biology.
With a single batch level the tool refuses — nothing to integrate.

## py_markers

`rank_genes_groups` (wilcoxon by default) on `layers.lognorm` when present, with `pts`
so the Seurat route's `min_pct` filter has an analogue. Tables are filtered by
`logfc_threshold` / `min_pct` / `only_positive` and sorted by adjusted p-value,
one CSV per cluster. The heatmap shows the top five genes per cluster.

## py_annotate

The same scoring rule as the Seurat route: per-cluster mean expression of the
panel genes minus the mean of all other genes, on the pre-scale snapshot.
Difference in means keeps panels of different sizes comparable. `min_score`
gates acceptance; below it a cluster stays `unassigned`. An explicit mapping
JSON overrides scoring but must label every cluster. The dot plot is the
evidence to check before trusting a label.

## py_pseudobulk

Sums the counts layer per sample (or sample × group), drops units under
`min_cells`, keeps genes expressed in at least `min_cells_expressing` units,
and writes the wrapper-compatible layout: `gene_id` first column, unit-named
columns, sample sheet whose `sample` column matches them. `wrapper_hint`
names exactly what to pass the DESeq2/limma/edgeR wrapper. Without a group
column there is no contrast to test and the hint is null.

## py_trajectory

PAGA coarse graph over clusters plus diffusion pseudotime (`dpt`) rooted at
`root_cluster`. Pseudotime is relative to the root and undefined (inf) for
cells in disconnected components; both facts are reported. Choosing the root
is a biological decision for the agent, not a default to accept silently.

## Lab mapping

| NBIS lab | Preset | Chain coverage |
| --- | --- | --- |
| 1 QC | `module: 1` | inspect → qc |
| 2 Dimensionality reduction | `module: 2` | + preprocess |
| 3 Integration | `module: 3` | + integrate, batch-compare |
| 4 Cluster markers | `module: 4` | + markers |
| 5 Sample-level DE | `module: 5` | + pseudobulk (wrapper handoff) |
| 6 Annotation | `module: 6` | + annotate, assign-labels |
| 7 Trajectory | `module: 7` | + trajectory |

The labs' optional sections that need packages this environment omits (scVI,
CellTypist models, Scanorama, gseapy, spatial) are reported by
`sc_py_env_probe`, never silently skipped.