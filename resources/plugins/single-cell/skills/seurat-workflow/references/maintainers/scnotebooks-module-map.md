# Module map

Per-module detail for the 13 scNotebooks modules, with the packages each step
needs and what exists in the plugin environment (`plugin:seurat`) today.
Verified against `environments/seurat/environment.yml` in this plugin.

## What the plugin environment provides directly

`Seurat` 6.0.0, `SeuratObject`, `SingleCellExperiment`, `ComplexHeatmap`,
`jsonlite`, `ggplot2`, `patchwork`, and the tidyverse (see the environment spec
for the full list). `Matrix` and `uwot` arrive as Seurat imports rather than as
explicit specs. `harmony` is **not** shipped.

Everything else in the tables below is **absent** from the plugin environment.
Run `sc_env_probe` to confirm on the machine in front of you rather than trusting
this table.

Since 0.5.0 the plugin also ships a Python/scanpy route (`scanpy-workflow`,
`plugin:scanpy`). Some R-route gaps have Python counterparts there — doublets
(Scrublet), batch integration (Harmony/ComBat), trajectory (PAGA + DPT) — but the
tables below stay R-environment facts: `sc_env_probe` computes them against
`plugin:seurat`, and the Python route reports its own gaps via `sc_py_env_probe`.

## Module by module

### 01 — Notebooks and databases (Python, 2 code cells)

Cell structure, Colab, and a tour of single-cell repositories: Single Cell
Expression Atlas, HCA Data Portal, CELLxGENE, PanglaoDB, CellTypist, GEO, SRA.

*In Phi:* the notebook concepts belong to Phi's own notebook support; the
databases are covered by the `transcriptomics-single-cell` connector
(`cellxgene_list_datasets`, `hca_search_projects`, `cellguide_cell_type`,
`gtex_gene_expression`). Nothing to port as a script.

### 02 — R and ggplot2 (R, 54 code cells)

R syntax, vectors and data frames, then the grammar of graphics.

*In Phi:* the `Visualization` agent with `omics-visualization` for the plotting
part. No port needed; keep the notebook as reading.

### 03 — Raw reads and Cell Ranger (bash, 29 cells)

Shell basics (`ls`, `cd`, `mv`, `cat`, `head`, `wget`), SRA toolkit, downloading
and running `cellranger count`, building a chr22 reference with `samtools faidx`
and `cellranger mkref`.

*In Phi:* belongs to the wrapper tree, not a skill. The repository already ships
148 nf-core modules plus the `rna-seq` workflow, and the `install-nextflow`,
`run-module`, `launch-workflow`, and `create-workflow` skills drive them. Cell
Ranger itself is a licensed 10x binary and is not redistributable; the notebooks
download it from 10x at run time, which a read-only environment cannot do.

### 04 — QC, clustering, annotation, functional analysis (R, 141 cells)

The largest module, and the one this plugin ports. Sections:

| Section | Needs | In plugin env |
| --- | --- | --- |
| Read 10x `.h5`, sparse matrix handling | `hdf5r` (for `.h5`) | no → read the mtx directory instead |
| QC metrics, violin/histogram/scatter | Seurat, ggplot2 | yes |
| `scDblFinder` doublets | `scDblFinder`, `SingleCellExperiment` | SCE yes, scDblFinder no |
| Normalization, variable genes, scaling, PCA, UMAP, clusters | Seurat | yes |
| `SCTransform` | Seurat + `glmGamPoi` (optional) | Seurat yes |
| Markers, `FindAllMarkers` | Seurat (+ `presto` for speed) | yes, slower |
| `SingleR` + `celldex` annotation | `SingleR`, `celldex` | no → panel scoring instead |
| Manual marker panel + dot plot | Seurat | yes |
| Milo differential abundance | `miloR` | no |
| GO enrichment | `clusterProfiler`, `org.Hs.eg.db` | no |
| SoupX ambient RNA | `SoupX`, `DropletUtils` | no |

### 05 — Multi-sample integration (R, 36 cells)

`ifnb` dataset, `SplitObject`, `FindIntegrationAnchors`, `IntegrateData`, then
Harmony, then `MixingMetric` / `LocalStruct` comparisons.

*In Phi:* the normalization/scale/PCA/cluster half runs through
`sc_preprocess`. Anchor integration and Harmony kernels need
`SeuratWrappers`/`harmony` handling that the plugin does not expose as a script
surface; a project environment could add it. `ifnb` comes from `SeuratData`,
which downloads on first use — another read-only-environment blocker.

### 06 — Trajectory and pseudotime (R, 24 cells)

Monocle3: `as.cell_data_set`, `cluster_cells`, `learn_graph`, `order_cells`,
`plot_cells`, `graph_test`.

*In Phi:* **blocked.** `monocle3` is not in the plugin environment. The notebook also
installs it from GitHub inside the session. Report the gap.

### 07 — Cell-cell communication (R, 29 cells)

LIANA: `liana_wrap` across cellphonedb/natmi/sca/logfc, `liana_aggregate`,
`liana_dotplot`, `heat_freq`. Installs `Seurat` 4.4.0 explicitly, which a
read-only environment cannot honour.

*In Phi:* **blocked.** `liana` is not in the plugin environment.

### 08 — Multimodal / CITE-seq (R, 11 cells)

CBMC dataset: RNA and ADT assays, CLR normalization, `FeatureScatter` between a
protein and its RNA, negative-control mouse cells.

*In Phi:* the workflow is ordinary Seurat, so `sc_preprocess` and
`sc_markers` apply once the object exists. The notebooks read `.h5` and
dense CSVs from GEO; a `.csv` count table works with `sc_inspect`, `.h5`
needs `hdf5r`.

### 09 — TCR sequencing (R, 18 cells)

10x V(D)J: `Read10X` of a filtered feature-barcode matrix, QC, clustering, then
merging the `*_contig_annotations.csv` tables and plotting clonotype frequency.

*In Phi:* the gene-expression half runs through the seurat chain. Clonotype
merging is a small table join the notebooks do inline; it is not part of this
port. CELLxGENE/HCA do not distribute V(D)J tables.

### 10 — Spatial transcriptomics (R, 59 cells)

Visium (10x breast cancer sections), `semla`: `ReadVisiumData`, `LoadImages`,
`MapFeatures`, NMF via `singlet`/`RcppML`, `RunNNLS` deconvolution, Enrichr.

*In Phi:* **blocked as written.** `semla` and `singlet` are not in the plugin environment, and the
section downloads ~2 GB from 10x per run. Spatial is genuinely out of scope for
the current environment; say so plainly.

### 11 — scATAC-seq (R, 67 cells)

ArchR end to end: `createArrowFiles`, doublets, `addIterativeLSI`, Harmony,
clusters, gene scores, peak calling with MACS2, `addMotifAnnotations`, chromVAR,
co-accessibility, footprints, trajectories.

*In Phi:* **blocked.** `ArchR` is not in the plugin environment, and MACS2 is a separate
Python/binary dependency the module installs with `pip`.

### 12 — Alternative polyadenylation (Python, 32 cells)

SCAPE-APA: conda install, `gen_utr_annotation` on a GFF3, `prepare_input` on a
BAM, `infer_pa`, `merge_pa`, plus pandas/matplotlib plots.

*In Phi:* out of scope for an R workflow. It is a CLI tool over BAM/GFF, closer
to the wrapper tree than to a skill, and the GFF it downloads is large.

### 13 — FAIR data sharing (8 cells)

Metadata for GEO/SRA/BioProject/BioSample, file formats, MAGE-TAB for
ArrayExpress/Single Cell Expression Atlas, HCA submission tiers.

*In Phi:* reference reading. No executable step, and none is needed.

## Score card for the current environment

Runnable core: 01 (reading), 02, 03 (via wrappers), 04 (core), 05 (partly),
08, 09, 12 (out of scope), 13 (reading).

Blocked by a missing package: 04 (Milo, SoupX, GO, scDblFinder), 06, 07, 10, 11.

That split is what `sc_env_probe` reports, computed rather than asserted.
