# single-cell (single-cell analysis) · 0.6.0

Two parallel single-cell RNA-seq toolchains: **R/Seurat 6** and **Python/scanpy 1.13**, each
with a private locked environment, routed by user preference. Includes a map of the 13 teaching
modules of [scNotebooks](https://github.com/integrativebioinformatics/scNotebooks); the Python
route's notebook presets adapt the scanpy labs of the
[NBIS workshop-scRNAseq](https://nbisweden.github.io/workshop-scRNAseq/home_contents.html).

Two parallel single-cell toolkits (R/Seurat 6 and Python/scanpy 1.13), each with its
own locked private environment, routed by user preference. Includes a map of the 13
scNotebooks teaching modules; Python notebook presets adapt the NBIS scanpy labs.

## Components

| Component | Role |
|--|--|
| Agent `SingleCell` | Role, boundaries, and the **R/Python route rules**; binds no environment and builds only the chosen route's environment |
| Skill `seurat-workflow` | The R route's executable chain: 8 `sc_*` tools (including notebook rendering) |
| Skill `scanpy-workflow` | The Python route's executable chain: 10 `sc_py_*` tools (including notebook rendering) |
| Reference `scnotebooks-curriculum` | Module-by-module map of the 13 modules, tucked into `seurat-workflow/references/maintainers/`; documentation only |
| Environment `seurat` | R 4.5 + Seurat 6.0.0 + IRkernel, locked for three platforms |
| Environment `scanpy` | Python 3.13 + scanpy 1.13 + ipykernel, locked for three platforms |

The two skills run in `plugin:seurat` and `plugin:scanpy` respectively and do not depend on a
shared base environment. Specs and the three-platform locks live in `environments/<name>/`.
The notebook kernels are named `plugin-single-cell-seurat` / `plugin-single-cell-scanpy` and
are synchronized by the Phi host once the environment is ready.

## Route selection (the agent's job, not the renderer's)

- The user explicitly says R/Seurat or Python/scanpy → take that route and never mix.
- Only files given: `.rds` → R route; `.h5ad` → Python route; only a 10x matrix directory → ask about their preference.
- Switching routes after one is finished is a new task; cross-route handoff goes only through
  two interfaces: the `pseudobulk` counts table (handed to the `deseq2/limma/edger` wrapper) or
  an explicit object-conversion request.
- The full rules are in `agents/SingleCell.md` and the first Rules item of both SKILL.md files.

## Tools

R route (prefix `sc`, environment `plugin:seurat`):

| Tool | Purpose |
|--|--|
| `sc_inspect` | Read a `.rds`, a 10x matrix directory, or a counts table, report shape, assay/layer, reductions, and metadata columns, and list the steps it can already do |
| `sc_qc` | Mitochondrial/ribosomal/hemoglobin fractions, doublet detection, threshold filtering; reports how many cells each filter removed and writes violin plots |
| `sc_preprocess` | Normalization (LogNormalize or SCTransform), variable features, scale, PCA, neighbor graph, multi-resolution clustering, UMAP/t-SNE |
| `sc_markers` | `FindAllMarkers`, one CSV per cluster plus a marker heatmap |
| `sc_annotate` | Marker-panel scoring annotation or an explicit mapping; writes a dot plot |
| `sc_pseudobulk` | Per-sample (or sample×group) aggregate counts table and sample sheet, with column names following the DE wrapper convention, and a returned `wrapper_hint` |
| `sc_env_probe` | Reports missing packages and the parts that cannot run, module by module |
| `sc_notebook` | Renders an R notebook: preset modules 4/5/8/9 or a steps array, dependencies filled in automatically, parameters gathered in a config cell, injectable/mergeable |

Python route (prefix `sc_py`, environment `plugin:scanpy`):

| Tool | Purpose |
|--|--|
| `sc_py_inspect` | Read a `.h5ad`, a 10x `.h5`/matrix directory, or a counts table, report shape, layers (including the lognorm snapshot), obsm, and obs columns |
| `sc_py_qc` | QC metrics, Scrublet doublets, threshold filtering; reports how many cells each filter removed and writes distribution plots |
| `sc_py_preprocess` | normalize+log1p (keeping a lognorm snapshot), seurat_v3 variable genes, scale, PCA, Leiden, UMAP/t-SNE |
| `sc_py_integrate` | Harmony or ComBat batch integration, before/after UMAP comparison, re-clustering after integration |
| `sc_py_markers` | `rank_genes_groups`, one CSV per cluster plus a marker heatmap |
| `sc_py_annotate` | Marker-panel scoring annotation or an explicit mapping; writes a dot plot |
| `sc_py_pseudobulk` | The sample-level counts table, isomorphic to the R route (expression filtering works on unit counts, same semantics), with a returned `wrapper_hint` |
| `sc_py_trajectory` | PAGA abstract graph plus diffusion pseudotime (DPT), with a chosen root cell |
| `sc_py_env_probe` | Reports missing packages in the scanpy ecosystem (bbknn, scVI, CellTypist, and so on) and the parts that cannot run |
| `sc_py_notebook` | Renders a Python notebook: preset NBIS labs 1–7 or a steps array; the remaining capabilities are the same as the R renderer |

Both routes' `build_notebook.py` is **the same byte-identical renderer** (asserted by the
repository tests), and the language is decided by the `--language` / `language` argument; one
kernel per notebook, and a language mismatch is an explicit error.

## Two modes (both available on each route)

**Inspect as you go (learning, exploration)**: render a notebook, then run it cell by cell with
Phi's `notebook_run_cell`.

```
sc_notebook(module=4, data="tenx_dir/", self_check=true)          # R
sc_py_notebook(module=3, data="batched.h5ad", self_check=true)    # Python
notebook_read → notebook_run_cell → the user edits the config cell → only the affected cells re-run
```

The template is source material to render, not an object to copy: the fragment library
(`assets/fragments/`) is assembled, using the steps and parameters the agent passes in, into a
notebook with its placeholders filled, and lands in the user's working directory. Which steps to
pick and which thresholds to use is the agent's judgement based on the user's question and the QC
distributions (the rules are in SKILL.md); the renderer only does deterministic filling. All
thresholds, paths, and column names are gathered in the config cell in the first cell; render
provenance (steps, parameter sources, fragment hashes, data fingerprint) is recorded in the
notebook metadata. `self_check` really runs every cell in the output directory, headless
(Agg + autoshow=False), and never pops up a native plot window.

**Run it through in one go (reproducible)**:

```
inspect → qc → preprocess → [integrate] → markers → annotate → pseudobulk
```

Each step reads one object and writes one (`.rds` / `.h5ad`), can be re-run on its own, and
records its parameters in the returned JSON. On the Python route, `trajectory` can follow
markers. The two routes and the two carriers should reach the same conclusions; if they do not,
one of them has a bug.

**Raw data processing is out of scope**: alignment and Cell Ranger belong to the wrapper tree.
Notebooks and script chains both start from a count matrix or an existing object.

## A single analysis source and a knowledge layer (0.6.0)

Each route has only **one copy of the analysis code**: `assets/fragments/` (flat, one language;
presets and parameter defaults in `assets/steps.json`). The notebook's cells are the fragment
text itself; the script tools are thin adapters — `scripts/_fragments.{R,py}` runs the same batch
of fragments in one global environment, the adapter only maps CLI flags to configuration globals
named like the config cell's, and `scripts/reports/*` only serializes JSON. The fragments'
teaching output goes to stderr, leaving a single JSON document on stdout (R uses
`sink(stderr())`, Python uses `redirect_stdout`). Checked field by field before and after the
refactor: the R 7/7 and Python 10/10 reports match the old contract (only the output paths
differ).

Knowledge sits where it is used: one `references/step-notes/<step>.md` per step (purpose / how it
works / reading the results / pitfalls — the notebook's explanation cells, paired with the
fragments by step id); each route has an
[analysis-playbook.md](skills/seurat-workflow/references/analysis-playbook.md) — the agent's
decision manual, organized per stage as "call → read (the exact JSON fields and figures) → decide
(starting parameters + adjustment rules) → red lines", with the theory compressed to one line and
sources and licences at the end (OSCA CC BY 4.0, sc-best-practices Apache 2.0, the Cambridge
course GPL-3, NBIS CC-BY-4.0). `references/README.md` indexes by function: agent decisions
(playbook, step-notes, cross-route handoff) are kept apart from maintainer provenance/comparison
material (`maintainers/`). All of it is original prose; no upstream text is copied.

## Relationship to the upstream courses

- **scNotebooks** (CC BY 4.0, Rojas-Hidalgo et al., *Nat Genet*, 2026): the course is not a
  pipeline (Colab installs packages, `/content/` paths are hard-coded, state stays in the kernel).
  The plugin rewrites the runnable parts as scripts and fragments and says plainly what cannot run.
  Runnable: the module 04 main line, half of the 05 preprocessing, the expression parts of 08/09;
  not runnable: 06 (monocle3), 07 (liana), 10 (semla), 11 (ArchR) — the R environment lacks those
  packages; 03 (Cell Ranger) belongs to the wrapper tree and 12 (SCAPE-APA) to separate tools. See
  the [module map](skills/seurat-workflow/references/maintainers/scnotebooks-module-map.md).
- **NBIS workshop-scRNAseq** (CC-BY-4.0, as stated in the repository LICENSE; the site footer's
  GPL-3 does not match the repository): the scanpy labs serve as a **reference** only, and their
  text is not redistributed. The fragment code is original to this plugin and semantically
  aligned with labs 1–7 (QC, dimensionality reduction, integration, clustering+markers,
  sample-level differential expression, annotation, trajectory); each preset's intro and
  SOURCE.json record the source and licence. See the
  [NBIS porting notes](skills/scanpy-workflow/references/maintainers/nbis-porting.md) for the
  porting decisions.

## Boundaries

- Upstream repository content is not modified, and their notebooks, images, or example data are
  not redistributed.
- No packages are installed into managed environments. A gap is reported as a gap: the R route
  lacks monocle3/liana/semla/ArchR/scDblFinder; the Python route lacks bbknn (no py313 build on
  bioconda), scVI/CellTypist/Scanorama/gseapy/spatial. `sc_env_probe` / `sc_py_env_probe` report
  them module by module.
- The Python route can do trajectory (PAGA+DPT) and doublets (Scrublet); the R route still marks
  the corresponding modules as missing packages and does not pretend to run them.
- Spatial transcriptomics and scATAC-seq are out of scope.

## Why not a Nextflow wrapper

Making the single-cell workflow into a wrapper is not impossible, but it does not fit, because the
execution model differs:

- **A wrapper runs with `-profile conda` and resolves dependencies at run time**; **the plugin
  environment is a lock plus read-only** (the client never resolves dependencies). Object
  structures are incompatible between Seurat/scanpy major versions (this plugin locks Seurat
  6.0.0 / scanpy 1.13), and handing that to a solver on another day or another machine may
  resolve different versions.
- **Mismatched granularity**: single-cell analysis is human-in-the-loop (set thresholds after
  seeing the distributions, annotate after seeing the markers), and parameters change between
  steps once results are read; a script tool call completes one step and returns a structured
  result.
- **The one step that should use a wrapper does**: both routes' `pseudobulk` turn the object into
  input the `deseq2/limma/edger` adapter can consume and return a `wrapper_hint`; the R route's
  interface has been verified against the real DESeq2 wrapper.

See sections 8 and 10 of the
[integration design](../../../docs/design/scnotebooks-integration.zh-CN.md) for the full
argument.