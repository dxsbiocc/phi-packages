# OSTA chapter mapping (maintainers)

[OSTA](https://bioconductor.org/books/OSTA/) (Orchestrating Spatial Transcriptomics Analysis
with Bioconductor, Crowell et al., bioRxiv 2025, doi:10.1101/2025.11.20.688607) is an R/Bioconductor
book on spatial analysis, 40 chapters in v1.2.2. This plugin takes the **Python route only**
(scanpy + squidpy + spatialdata) and translates OSTA's analysis steps and principles into the
Python ecosystem; R package names serve only as conceptual anchors, not as dependencies. Licence
CC BY 4.0 (pending final check of DESCRIPTION). The knowledge is an original distillation; no
upstream text/figures/data are copied.

## Background (1–8)

| OSTA chapter | R anchor | Implementation in this plugin |
| --- | --- | --- |
| 4 Infrastructure | SpatialExperiment / SingleCellExperiment | AnnData (`obsm["spatial"]` coordinates, `uns["spatial"]` image); spatialdata `SpatialData` for multi-sample/multi-scale |
| 6 Importing data | TENxVisiumData, SpatialExperiment reading | The `load` fragment + `spatialdata-io` (visium/visium_hd); a custom reader for Stereo-seq GEM/GEF |
| 8 Python interoperability | reticulate calling Python | This plugin is itself the Python route, so no interoperability layer is needed |

## Sequencing-based platforms (9–17)

| OSTA chapter | R anchor | Implementation in this plugin |
| --- | --- | --- |
| 10 Reads to counts | Cell Ranger / spaceranger | The wrapper tree (not in the notebook) |
| 11 Quality control | per-spot QC, tissue coverage | `qc-metrics`/`qc-plots`/`filter` (distributions on spatial coordinates and tissue maps) |
| 12 Intermediate processing | normalize/HVG/PCA/cluster | `normalize`/`hvg`/`scale`/`pca`/`cluster` |
| 13 Deconvolution | RCTD(spacexr), cell2location, SPOTlight, stereoscope, Tangram | `deconvolve` (scipy NNLS, needs an scRNA reference; cell2location/tangram are project-environment gaps) |
| 14–15 Visium DLPFC / CRC workflow | full workflow | presets (the standard Visium chain) |
| 16–17 Visium HD bin / seg | square bin / cell segmentation | `load` recognises HD bin and seg output; the analysis units are the bin/cell respectively |

## Imaging-based platforms (18–26) — not covered in v0.1

Xenium/CosMx/MERFISH are not in v0.1 (Visium + Stereo-seq only). Their concepts of
**neighbourhood analysis (22)** and **cell communication (23)** are reused by the
sequencing-based platforms: `nhood`, `communication`. Segmentation (19) corresponds to
Stereo-seq cellbin — segmentation is done upstream by SAW and the plugin only reads the output.
Sub-cellular analysis (24) and image deep learning (33–34) are gaps.

## Platform-agnostic analysis (27–34)

| OSTA chapter | R anchor | Implementation in this plugin |
| --- | --- | --- |
| 27 Normalization | logNormCounts, sctransform | `normalize` (library-size + log1p) |
| 28 Dimensionality reduction | PCA, UMAP | `pca`/`cluster` (UMAP) |
| 29 Clustering & annotation | Leiden, marker panel | `cluster`/`markers`/`annotate`/`assign-labels` |
| 30 Feature selection & testing | nnSVG, SPARK, SpatialDE (SVG) | `svg` (squidpy `spatial_autocorr`: Moran's I / Geary's C) |
| 31 Feature sets & signatures | panel scoring, GSVA | `annotate` (panel-mean-minus-rest) |
| 32 Spatial statistics | neighbourhood enrichment, co-occurrence, Moran | `nhood` (squidpy `nhood_enrichment`/`co_occurrence`) |
| 33–34 Image analysis / deep learning | image segmentation/deep learning | gap (project environment) |

## Multi-sample analysis (35–37)

| OSTA chapter | R anchor | Implementation in this plugin |
| --- | --- | --- |
| 35 Differential spatial patterns | cross-condition SVG differences | `spatial-diff` (differences in spatial patterns between conditions) |
| 36 Differential colocalization | cross-condition colocalisation differences | `spatial-diff` (colocalisation comparison) |
| 37 Structure-based analysis | tissue structure comparison | partly covered by `nhood`; the rest is a gap |

## Cross-platform analysis (38–40) — gaps

Registration (38), Imputation (39) and Xenium×Visium (40) are not in v0.1 and are listed as
future work/project environments.

## Division of labour with the single-cell plugin

The QC/normalize/cluster/marker/pseudobulk semantics of OSCA (single-cell) are already implemented
in the single-cell plugin; this plugin reuses its threshold semantics and panel scoring rules. What
is new is the **spatial dimension**: coordinates, adjacency graphs, SVG, neighbourhoods,
deconvolution, communication and cross-sample spatial patterns. The scRNA reference for
deconvolution comes from the single-cell plugin.