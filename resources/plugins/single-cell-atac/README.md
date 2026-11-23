# single-cell-atac

Single-cell chromatin accessibility (scATAC) analysis plugin. **Python single route**: snapatac2 + scanpy, private locked environment
`plugin:scatac`, tool prefix `scatac_`. The main engine is **snapatac2** (Rust backend, file-backed AnnData,
plotly figures); scanpy is used only downstream of the gene activity matrix. The analysis unit is the
**cell (barcode)**, not a spot/bin.

| Input starting point | Source | Notes |
| --- | --- | --- |
| `fragments.tsv.gz` | 10x Cell Ranger ATAC | Preferred; snapatac2 imports it directly and computes QC/TSSe |
| BAM | alignment output | Convert to fragments with `make_fragment_file` first, then import |
| `.h5ad` | existing snapatac2/scanpy object | backed or in-memory |
| peak×cell matrix | csv/tsv/h5/10x | features already present; skips the tile matrix |

**Raw read alignment / Cell Ranger ATAC is not done here**; that is the job of the nf-core wrapper.

## Design principles

The same contract as the single-cell / spatial plugins:

- **The plugin only does deterministic filling/rendering.** Intent parsing, step selection and threshold decisions
  belong to the agent (the rules live in
  [`agents/SingleCellAtac.md`](agents/SingleCellAtac.md) and the skill's SKILL.md).
- **The fragment library is the single analysis source.** The teaching code in
  `skills/scatac-workflow/assets/fragments/` drives both the script tools (thin adapters + a JSON report layer)
  and the notebook renderer; the two forms share one copy of the code.
- **Private locked environment; do not install packages inside it.** Missing packages or genome references are
  reported as facts by `scatac_env_probe`; the next step is a project environment (env_request) or instructions
  to download into the project directory.
- **The notebook is real code cells + per-step Markdown teaching**, rendered into the user's working directory,
  not a copied template directory.

## Two hard constraints of snapatac2

- **backed AnnData**: the object is backed by the `.h5ad` file and syncs automatically; **`close()` is mandatory
  before the process exits**, otherwise the HDF5 is corrupted. The script chain writes and closes at every step;
  the notebook ends with a close cell.
- **plotly figures, not matplotlib**: `snap.pl.*` returns plotly Figures, and the script chain saves them as
  self-contained `.html`. matplotlib is used only when the gene activity downstream goes through scanpy.

## Analysis chain

`scatac_inspect` → `scatac_qc` (fragment count/duplication rate/mitochondrial fraction/TSSe + threshold distributions) → `scatac_preprocess`
(tile matrix/feature selection/scrublet doublets/spectral embedding/Leiden clustering/UMAP) → `scatac_peaks`
(per-cluster MACS3 peaks + merged peak set + peak-count matrix) → `scatac_markers` (per-cluster
marker regions, optional pairwise differential accessibility `run_diff`) → `scatac_annotate` (gene activity +
panel/scRNA reference/explicit mapping labels) → `scatac_motif` (motif enrichment on marker regions) →
`scatac_integrate` (multi-sample batch integration, harmony/scanorama/mnc on the spectral embedding).
`scatac_notebook` is the parallel teaching path; `scatac_env_probe` probes before you promise heavy-dependency capabilities.

## Genome references

- **chrom_sizes bundled**: snapatac2 ships hg38/mm10/hg19/mm39 and works offline.
- **blacklist inlined**: the ENCODE hg38/mm10 blacklist is inlined as base64-gzip in
  `assets/fragments/blacklist.py`; it is decoded offline and needs no external file.
- **Gene annotation (GFF3 ~50MB) and genome FASTA (~845MB) are lazy downloads**, not distributed with the environment: the user supplies
  `input-path`, or downloads them once into the project directory; when missing, report the gap — never
  download large references into the read-only plugin environment.
- **The motif library** (JASPAR/cis_bp) is a snapatac2 download resource, not distributed with the plugin.

## Knowledge layer

`skills/scatac-workflow/references/` is organized by function:

- [`analysis-playbook.md`](skills/scatac-workflow/references/analysis-playbook.md) — the agent's
  operating manual: which tool to call at each step, which JSON fields and figures to read, how to set
  thresholds, and what the red lines are.
- [`red-line-cheatsheet.md`](skills/scatac-workflow/references/red-line-cheatsheet.md) —
  12 hard red lines × common misreadings quick-reference table, for use at runtime.
- [`step-notes/`](skills/scatac-workflow/references/step-notes/) — the purpose/principle/interpretation/pitfalls of each step,
  i.e. the notebook's teaching cells.
- [`maintainers/implementation-notes.md`](skills/scatac-workflow/references/maintainers/implementation-notes.md)
  — implementation notes and provenance (fragment library architecture, backed AnnData constraints, snapatac2 API mapping, JSON contract).

The knowledge is an original distillation of the [SnapATAC2](https://github.com/kaizhang/SnapATAC2) (MIT) documentation and
[sc-best-practices](https://www.sc-best-practices.org/) (Apache 2.0); the ENCODE blacklist comes from
doi:10.1038/s41598-019-45839-z. No upstream text/figures/data are copied.

## Environment

`environments/scatac/environment.yml` declares the dependencies (pure conda, no pip), and `locks/` holds the lock files for the three
platforms (darwin-arm64, darwin-x64, linux-x64). Solving goes through the tuna mirror, then the host is rewritten back to the official one.
macs3 is distributed with the environment as a snapatac2 dependency.

## Gaps (reported as facts by `scatac_env_probe`, handled via a project environment)

- **chromVAR/ArchR/Signac**: the R/Bioconductor route, not implemented. The plugin's TF layer is motif enrichment (sequence prediction),
  **not chromVAR per-cell TF deviations**, and motif hits must not be treated as measured TF binding either.
- **episcanpy**: although it is on bioconda, at runtime `from numpy import NaN` is incompatible with numpy 2, so it is not distributed with the environment.
- **scvi-tools / magic-impute**: pip-only, and the locked environment does not support pip.
- **TF-gene GRN**: not implemented (a snapatac2 inject extension, not part of this plugin).

## Relationship to other plugins

- **single-cell**: the scRNA reference for label transfer / gene-activity annotation comparison comes from the SingleCell agent,
  via the shared exit (counts+obs exchange / pseudobulk table). Gene activity is a **proxy** for accessibility, not measured RNA.
- **wrappers**: raw reads to fragments (Cell Ranger ATAC, alignment) is the nf-core wrapper's job, not a notebook cell.
- **chromatin-profiling** (planned sibling plugin): bulk ChIP/CUT&Tag is another plugin (R route, purely downstream,
  starting from called peaks/counts), kept separate from this plugin by analysis paradigm.