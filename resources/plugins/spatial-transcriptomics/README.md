# spatial-transcriptomics

Spatial transcriptomics analysis plugin. **Python route only**: scanpy + squidpy +
spatialdata in the private locked environment `plugin:spatial`, tool prefix `st_py_`.
Covers two sequencing-based platforms:

| Platform | Input | Analysis unit |
| --- | --- | --- |
| 10x Visium | Cell Ranger matrix directory/.h5 + spatial coordinates + tissue image | spot (holds several cells) |
| 10x Visium HD | square bin (8/16µm) or seg (cell segmentation) output | bin or cell |
| BGI Stereo-seq | GEM/GEF (SAW output) | square bin (binN) or cellbin (cells) |

## Design principles

Same contract as the single-cell plugin:

- **The plugin only does deterministic filling/rendering.** Intent parsing, step
  selection and threshold decisions belong to the agent (the rules live in
  [`agents/SpatialTranscriptomics.md`](agents/SpatialTranscriptomics.md) and the skill's SKILL.md).
- **The fragment library is the single source of analysis code.** The teaching code in
  `skills/spatial-workflow/assets/fragments/` drives both the script tools (thin adapters
  + JSON report layer) and the notebook renderer; the two forms share one copy of the code.
- **Private locked environment; do not install packages in it.** A missing package is
  reported as a fact by `st_py_env_probe`, and the next step is a project environment.
- **A notebook is real code cells + per-step Markdown teaching**, rendered into the user's
  working directory, not a copied template directory.

## Analysis chain

`st_py_inspect` → `st_py_qc` → `st_py_preprocess` → `st_py_markers` →
`st_py_annotate` → `st_py_spatial` (SVGs, neighbourhood enrichment/co-occurrence,
interaction matrix, centrality, Ripley's L, cross-sample differential spatial
patterns) → `st_py_niche` (UTAG-style spatial niche segmentation & cross-stage
dynamics) → `st_py_trajectory` (Scanpy PAGA + spatial diffusion pseudotime &
gene-program scoring) → `st_py_deconvolve` (needs an scRNA reference) →
`st_py_communication` (cell communication, optionally stratified by niche).
`st_py_notebook` is the parallel teaching path.

## Knowledge layer

`skills/spatial-workflow/references/` is organised by function:

- [`analysis-playbook.md`](skills/spatial-workflow/references/analysis-playbook.md) — the agent's
  operating manual: which tool to call at each step, which JSON fields and figures to read, how
  to set thresholds, and what the red lines are.
- [`step-notes/`](skills/spatial-workflow/references/step-notes/) — purpose/principle/interpretation/pitfalls
  for each step, i.e. the notebook's teaching cells.
- [`maintainers/`](skills/spatial-workflow/references/maintainers/) — implementation notes and
  provenance (OSTA chapter correspondence, Stereo-seq methodology).

The knowledge is an original distillation of [OSTA](https://bioconductor.org/books/OSTA/)
(Bioconductor, CC BY 4.0) and the spatial chapters of
[sc-best-practices](https://www.sc-best-practices.org/); the methodology follows 7 high-impact
Stereo-seq application papers. No upstream text/figures/data are copied.

## Environment

`environments/spatial/environment.yml` declares the dependencies, and `locks/` holds the lock
files for the three platforms (darwin-arm64, darwin-x64, linux-x64). Solving goes through the
tuna mirror, then the URLs are rewritten back to the official host.

## Relation to other plugins

- **single-cell**: the scRNA reference for deconvolution / label transfer comes from the
  SingleCell agent, through the shared gateway (counts+obs exchange / pseudobulk table).
  SingleCell lists spatial as a gap; this plugin fills it.
- **wrappers**: raw reads to counts (Cell Ranger, SAW, alignment) is the wrappers' job, not a
  notebook cell.