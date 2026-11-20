# scNotebooks curriculum, ported

[scNotebooks](https://github.com/integrativebioinformatics/scNotebooks) is a
single-cell and spatial omics teaching curriculum: 13 Jupyter notebooks in
English, Spanish, Portuguese, and Chinese, run on Google Colab or in a Docker
image. This reference is the index that maps it onto Phi; the executable work
is this skill's tool chain and fragment library.

## What the notebooks are, and what they are not

Each notebook is a **lesson**, not a pipeline. The pattern throughout:

- a setup cell that installs packages at run time (`install.packages`,
  `BiocManager::install`, `remotes::install_github`) and often configures
  apt/bspm through a downloaded shell script;
- teaching markdown that explains a concept, then one code cell per idea;
- paths such as `/content/` that assume the Colab runtime;
- questions left for the reader, with the answer in the next cell.

Phi environments are **read-only and prebuilt**, so the setup cell is the first
thing that cannot survive the move. Neither can a step that reads
`/content/`. Neither should: an install cell inside a managed environment is a
reproducibility hazard, not a feature.

## The 13 modules

| # | Module | What it teaches | In Phi |
| --- | --- | --- | --- |
| 01 | Notebooks and databases | Jupyter basics; CELLxGENE, HCA, GEO, SRA, PanglaoDB, CellTypist | Notebooks: Phi's own. Databases: the `transcriptomics-single-cell` connector |
| 02 | R and ggplot2 | R syntax, data frames, grammar of graphics | The `Visualization` agent and `omics-visualization` |
| 03 | Raw reads and Cell Ranger | Linux commands, SRA toolkit, `cellranger count` | Wrapper tree (`nf-core` modules); `install-nextflow`, `run-module` |
| 04 | QC, clustering, annotation, function | Seurat: doublets, QC, normalization, PCA, clusters, markers, SingleR, Milo, GO, SoupX | `seurat-workflow` core; Milo/SoupX/GO need packages the plugin environment omits |
| 05 | Multi-sample integration | Seurat anchors vs Harmony, mixing metrics | Core with `seurat-workflow`; Harmony via a project environment |
| 06 | Trajectory and pseudotime | Monocle3 graph learning, pseudotime | Blocked on R: `monocle3` is not in the plugin environment; the scanpy route offers PAGA + DPT |
| 07 | Cell-cell communication | LIANA, multiple methods, consensus | Blocked: `liana` is not in the plugin environment |
| 08 | Multimodal (CITE-seq) | RNA + ADT assays, CLR normalization | Core with `seurat-workflow`; needs `.h5` reading (`hdf5r`) or a 10x dir |
| 09 | TCR sequencing | 10x V(D)J contigs, clonotypes | Core with `seurat-workflow`; clonotype calls come from Cell Ranger |
| 10 | Spatial transcriptomics | Visium, semla, NMF, deconvolution | Blocked as written: `semla` is not in the plugin environment |
| 11 | scATAC-seq | ArchR, LSI, peak calling, motifs | Blocked: `ArchR` is not in the plugin environment |
| 12 | Alternative polyadenylation | SCAPE-APA on BAM + GFF | Out of scope: a Python/CLI tool, not an R workflow |
| 13 | FAIR data sharing | GEO/SRA/ArrayExpress/HCA submission, MAGE-TAB | Reference reading only |

"Blocked" means **not runnable in the `plugin:seurat` environment as it stands**, not
impossible forever: those steps need a project environment (`env_request`) that
adds the package. Say which it is.

## How to use this with a learner

1. Ask which module, and whether they have their own data or want the notebook's
   example data.
2. Run `sc_env_probe` for that module and state what will run.
3. For a runnable module, follow this skill: inspect, QC, preprocess,
   markers, annotate, pseudobulk. Keep the notebook's narrative: explain the
   concept, then run the step.
4. For a blocked module, give the notebook's explanation as reading, name the
   missing package, and say what a project environment would need.

## What to tell people about the notebooks themselves

They are CC BY 4.0 teaching material by the scNotebooks authors, cited as
Rojas-Hidalgo et al., *Nat Genet* (2026). When you use their text or figures,
keep the attribution. The notebooks remain the best place to *read* the
explanations; Phi is where the analysis becomes executable and reproducible.

See [scnotebooks-module-map.md](scnotebooks-module-map.md) for the per-module
detail, including which package each step needs.