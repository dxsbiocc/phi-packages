# Porting decisions

Why this plugin looks the way it does, and what was deliberately not copied from
scNotebooks.

## The notebooks are lessons, not pipelines

scNotebooks is a teaching curriculum: 13 notebooks, roughly 510 code cells in
English, mirrored in Spanish, Portuguese, and Chinese. The code cells are written
to be *read one at a time*, with the explanation in the markdown beside them.
Their code is byte-identical across the four language folders (verified: all 13
modules match), so a translation is a text change only.

Three properties follow, and each one rules out a mechanical conversion:

1. **Every notebook installs its own packages at run time.** 45 apt/bspm setup
   cells and 46 install calls across the corpus. Phi environments are immutable
   and prebuilt; an install inside one is a reproducibility hazard. Install
   cells are dropped, and a missing package becomes a reported capability gap.
2. **Paths assume Colab.** `/content/` appears in 7 cells, plus 38 `shell_call`
   helpers and 53 `!`-prefixed shell cells. Everything is re-expressed against a
   project directory the engine resolves.
3. **State lives in the kernel.** Each cell mutates a variable in memory
   (`sampleA`, `sc_merged`, `SpatialData`). A tool call has no kernel, so each
   step reads an `.rds` and writes one. That is also what makes steps re-runnable
   and auditable.

## Why a plugin and not 13 skills

A skill is one capability with its own scripts; a plugin is a composite of
agents, skills, and their environment. The scNotebooks content shares one
environment (`plugin:seurat`), one agent role (single-cell analysis), and one chain of
steps. Splitting it into 13 standalone skills would duplicate the object I/O and
the argument handling 13 times, and would put 13 near-identical entries in the
catalog.

The executable surface stayed one skill per language (`seurat-workflow`,
`scanpy-workflow`). The teaching map first shipped as a separate
`scnotebooks-curriculum` skill; since 0.6.0 it lives in
`scnotebooks-curriculum.md` and `scnotebooks-module-map.md` (same directory),
so a learner can ask "what does module 6 teach" without a second skill in the
catalog.

## What was deliberately not ported

| Not ported | Why |
| --- | --- |
| Colab install and apt setup cells | Read-only environment; the environment is the answer to "is this installed" |
| `/content/` paths, `shell_call` helpers | Tied to the Colab runtime |
| `.ipynb` files themselves | They are 1.1 MB of code plus ~30 MB of base64 images; the notebooks stay the reading material at their own URL |
| The 146 MB `scNotebooks-Resources/` folder | Git LFS pointers and a 147 MB tarball; sample data is fetched by the notebook, and this repo does not carry datasets |
| Module 12's SCAPE-APA CLI | A Python/BAM/GFF pipeline, not an R workflow; closer to the wrapper tree |
| Module 3's Cell Ranger invocation | Licensed 10x binary, not redistributable; belongs to the wrapper tree |
| The Spanish and Portuguese folders | Same code as English; the Chinese folder is the one with a maintained translation to reference |

## Environment reality, stated instead of hidden

The plugin maintains its own environment (`plugin:seurat`, R 4.5 + Seurat 6.0.0)
rather than borrowing the shared `phi-r`, so the script chain and the notebook
kernel stay on one pinned version. The single-cell stack it ships is deliberately
small. Checking the spec directly:

- present: `Seurat`, `SeuratObject`, `SingleCellExperiment`, `ComplexHeatmap`,
  `jsonlite`, `ggplot2`, `patchwork`;
- absent: `monocle3`, `liana`, `SoupX`, `DropletUtils`, `semla`, `ArchR`,
  `scDblFinder`, `pheatmap`, `clusterProfiler`, `org.Hs.eg.db`, `SingleR`,
  `celldex`, `miloR`, `clustree`, `presto`, `hdf5r`.

Conda availability was checked for the gaps, because "add it to the environment" is only an
option if a package exists: `r-monocle3` (bioconda, 0.1.3/0.2.0), `r-archr`
(bioconda, osx-arm64 only, 1.0.3), `r-soupx`, `r-clustree`, `r-harmony`,
`bioconductor-milor`, `r-presto`, `r-hdf5r` exist; `r-liana`, `r-semla`,
`r-singlet`, `r-seuratwrappers`, `bioconductor-dropletutils`, `r-cellchat` do
not (they would need `sourcePackages` with pinned commits).

So the port does three things instead of pretending:

1. **Implemented the core that works:** inspect, QC, preprocess/cluster, markers,
   annotate, pseudobulk — the whole Module 4 main line, Module 5's preprocessing
   half, and Modules 8 and 9's expression half.
2. **Made the gaps executable knowledge:** `sc_env_probe` computes, per
   module, which required packages are missing and which optional sections are
   unavailable. It is a tool call, not a claim in a document.
3. **Kept the curriculum honest:** the module map marks each module runnable,
   partly runnable, blocked, or out of scope, and names the package responsible.

## Behavioral changes from the notebooks

| Notebook behaviour | Here | Why |
| --- | --- | --- |
| Thresholds hard-coded after looking at plots | Arguments, with the distribution reported | A pipeline run must be reproducible without a human reading a histogram |
| `scDblFinder` always runs | `--detect-doublets`, off by default, best effort | The package is absent from the plugin environment; silently skipping would be worse than asking |
| Panel in a code cell | JSON panel file, difference-of-means score | Comparable across panels of different sizes; reviewable as data |
| Manual subsetting to drop cells | Per-filter removal counts and a kept fraction | An audit trail, plus a warning when a threshold is too strict |
| `FindAllMarkers` result printed | Tables per cluster, stored on the object | Reusable by the next step without recomputation |
| One sample per analysis | Pseudobulk per sample × group with a design report | Connects to the repository's differential-expression wrappers |

## Attribution

The scNotebooks notebooks are CC BY 4.0 (Rojas-Hidalgo et al., *Nat Genet*,
2026). This port re-implements their analysis steps as scripts; it copies no
notebook text, images, or data. Where an explanation is worth keeping, the
curriculum references point at the upstream notebook rather than reproducing it.
