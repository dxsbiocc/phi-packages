# Implementation notes — scATAC Python route

What each tool computes and how the layers divide the work. The chain is
deterministic: every step reads one object, writes one, and reports the
parameters it used. The analysis lives in `assets/fragments/` (per-step teaching
notes in `references/step-notes/`) — the tools are thin adapters over the same
fragments the notebook renders, so **adapt the fragment, not the script**. The
engine is snapatac2 2.10 (file-backed AnnData, Rust core); scanpy 1.12 is used
only for gene-activity downstream (normalize/log1p/pca/concat) and matplotlib
panels. MACS3 3.1 ships in the environment for peak calling.

## Architecture

```
assets/steps.json            registry: parameters (name/type/default/comment),
                             steps (id/title/requires), presets, kernel, footer
assets/fragments/<step>.py   self-contained teaching code, exec'd in a shared
                             namespace; the SINGLE analysis source for both paths
references/step-notes/<step>.md  the notebook's per-step Markdown teaching cell
scripts/_common.py           flag parsing, emit/fail/info, obs_col/obs_has/var_col,
                             detect_kind, sorted_levels, parse_flags
scripts/_config.py           config_namespace(): builds EVERY parameter global
                             from steps.json defaults + flag overrides; structural()
                             adds data_path/output_object/figures_dir/emit/tool_input
scripts/_fragments.py        run_fragments(steps, ns) execs fragments with stdout
                             redirected to stderr; run_report(name, ns) execs reports/
scripts/<tool>.py            thin adapter: parse flags → build namespace →
                             run_fragments → run_report; stdout = exactly one JSON
scripts/reports/<name>.py    serializes fragment bookkeeping into the JSON contract
scripts/env_probe.py         package/module probe (read-only env facts)
schemas/env-probe.json       the one declared tool output schema (env_probe only)
build_notebook.py            renderer (RENDERER_VERSION=3): reads fragments +
                             step-notes, emits config cell + code/Markdown cells
```

**Why `_config.config_namespace`.** A fragment references ~48 parameter globals
(`min_counts`, `bin_size`, `peak_qvalue`, `run_diff`, …). Building the namespace
from the registry guarantees every global exists with its declared default,
coerced to the declared type (`int`→flag_int, `number`→flag_num, `bool`→flag_bool,
`list`→flag_nums, else flag_str), then overridden by flags (flag name ==
parameter name) and structural values. This keeps the adapter namespace a
superset of the notebook config namespace — the two paths agree by construction,
and a fragment never hits a `NameError` because a flag was omitted.

**stdout discipline.** `run_fragments` redirects fragment stdout (print, scanpy
logs) to stderr; the report's `emit` writes exactly one JSON object to stdout. A
tool's stdout is always parseable JSON; progress is on stderr (`[scatac] …`).

**Bookkeeping.** `setup.py` initializes the shared globals; fragments leave
`figure_paths`, `analysis_steps`, `analysis_notes`, `qc_scored`, `qc_filters`,
`obj_before`, `adata`, `peak_adata`, `gene_adata`, `loaded_kind`, `genome_info`,
`cluster_columns`, `embeddings`, `peak_tables`, `merged_peaks`, `markers`,
`marker_tables`, `diff_results`, `annotation_panel`, `annotation_scores`,
`annotation_assignments`, `annotation_automated`, `motif_results`,
`integration_info` in the namespace. Reports read these; they never recompute
analysis.

## backed AnnData constraints (snapatac2-specific)

snapatac2 keeps a **file-backed** AnnData (Rust). These validated behaviors are
the reason for the helper functions and the separate-object design — violating
them panics the Rust core or silently corrupts the file:

- **`.obs` is lazy.** Supports `obs["col"].to_numpy()` and `"col" in adata.obs`,
  but NOT `.columns`, `.keys()` or iteration (those raise `PanicException`).
  `obs_has`/`obs_col` in `setup.py` (mirrored in `_common.py`) are the only safe
  access. obs column **assignment** works (`adata.obs[label] = array`).
- **`.uns` is `builtins.PyElemCollection`.** NO `.get()` method. Supports
  `"key" in adata.uns`, item access `adata.uns["key"]` (raises `RuntimeError`
  "No such key" if missing) and `.keys()`. `peaks.py` reads
  `adata.uns["macs3"] if "macs3" in adata.uns else {}` — never `.uns.get(...)`.
- **`.var` membership is unreliable.** `var_col` uses try/except item access
  rather than `"col" in adata.var`.
- **`to_memory()`** returns a standard `anndata.AnnData` with a pandas DataFrame
  obs — used in `annotate.py` (`_gene_memory`) for scanpy downstream.
- **`write(filename, backend="hdf5")`** exports a copy; **`.close()`** closes the
  backing; **`.isbacked`** is a PROPERTY, **`.is_closed()`** is a METHOD. Every
  backed object (`adata`, `peak_adata`, `gene_adata`) is closed before exit
  (`save.py` closes the peak/gene objects; the adapter/notebook closes `adata`).
- **The tile matrix cannot be resized in place.** `make_peak_matrix(...,
  inplace=True)` on a backed tile X FAILS; `inplace=False, file=path` WORKS and
  produces a separate backed peak object with obs (incl. leiden) copied. Same for
  `make_gene_matrix`. This is why peaks and gene activity are separate files.
- **`filter_cells` uses truthy bounds.** `min_tsse=None` OR `0` disables the gate
  (`if min_tsse:`); it accesses `obs["tsse"]` only when truthy. `filter.py`
  pre-checks the tsse column and nulls the bound when absent, so a missing
  annotation skips the gate with a note instead of panicking.

## Input loading

`detect_kind` (in `setup.py`, mirrored in `_common.py`) classifies the input;
`load.py` dispatches. Import always works on a **temporary backed file**, never
the user's input (a saved `.h5ad` is copied first so read-write backing cannot
mutate the original):

| kind | loader | notes |
| --- | --- | --- |
| fragment / fragment-like-table | `snap.pp.import_fragments` | chrom_sizes from `resolve_chrom_sizes`; `min_num_fragments` drops empty barcodes; populates `n_fragment`/`frac_dup`/`frac_mito` and `uns["reference_sequences"]` |
| bam | `snap.pp.make_fragment_file` → `import_fragments` | Tn5 shift + barcode tag first; records `fragments_from_bam` |
| anndata (.h5ad) | `shutil.copyfile` → `snap.read(backed="r+")` | resumes a saved object; genome_info source=saved-object |
| 10x-matrix-dir / 10x-h5 | `snap.read_10x_mtx` | features already present; **no fragments** → peaks/scrublet gap |
| peak-file / unknown | raises ValueError | raw FASTQ alignment belongs to a Cell Ranger ATAC wrapper |

`resolve_chrom_sizes` returns a bundled dict (`snap.genome.<hg38/mm10/hg19/
mm39>.chrom_sizes`, offline) or parses a user `chrom_sizes` file; an unsupported
genome raises rather than guessing. `resolve_blacklist` returns a user BED path,
else decodes the inlined base64-gzip ENCODE list (hg38 636 regions / mm10 3435
regions) to a temp `.bed.gz`, else None. The blacklist is **inlined**, not a
shipped asset file, because the rendered notebook runs in the user's project
directory with no plugin asset tree. `internal.intersect_bed` accepts `.bed.gz`
and plain `.bed`.

## inspect

Runs `setup → load → qc-metrics` (+ `save` when `output`). Reports `kind`,
`genome` provenance, `shape` (fragments import with 0 features), `metadata`
columns worth using (numeric min/max/median or category levels), and
`candidate_steps` computed from what the object already contains (`n_fragment` →
qc/preprocess/gene-activity; `X_spectral` → peaks/markers; `leiden` →
markers/annotate). `min_num_fragments` stays at its import default — per-cell
gating belongs to qc. `import numpy as np` is at the adapter top (the metadata
summary uses it outside the fragment namespace).

## qc

`qc_filter` runs `setup → load → qc-metrics → qc-plots` (plot_tag="before"),
then sets `plot_tag="after"` and runs `filter → qc-plots`, then `save`.
`qc-metrics.py` reads the import metrics, computes `snap.metrics.frag_size_distr`
(offline nucleosome ladder) and `snap.metrics.tsse` (needs `gene_anno`; skipped
with a gap note when null rather than a slow ~50 MB download), and snapshots
`obj_before` as plain numbers + per-sample counts BEFORE filtering (snapatac2
filters in place, so a reference would shrink too). `filter.py` applies
`snap.pp.filter_cells` with per-gate removed counts; scrublet is NOT here (it
needs the feature matrix, so it lives in `features.py`). The report emits
`distribution[]` (min/q01/median/q99/max per metric), `filters[]`,
`before`/`after` (`kept_fraction`), `features_scored`, `figures[]`, `warnings`.

## preprocess

`preprocess_cluster` runs `setup → blacklist → load → qc-metrics → qc-plots →
filter → features → embed → cluster → save`.
- `features.py`: `snap.pp.add_tile_matrix(bin_size)` when `n_vars==0` (a 10x
  matrix already has features); `snap.pp.select_features(n_features, blacklist)`;
  then `snap.pp.scrublet(features="selected", expected_doublet_rate)` +
  `snap.pp.filter_doublets(probability_threshold)` when `doublet_filter`.
- `embed.py`: `snap.tl.spectral(n_comps, features="selected")` →
  `obsm["X_spectral"]`; `snap.pp.knn(n_neighbors, use_rep="X_spectral")`. No
  variance-ratio curve (unlike PCA).
- `cluster.py`: `snap.tl.leiden(resolution)` per requested resolution (one →
  `obs["leiden"]`, several → `obs["leiden_<r>"]`); `snap.tl.umap(use_rep=
  "X_spectral")` + `snap.pl.umap(show=False)` when `run_umap`. The graph and
  Leiden never use UMAP.

Reports `steps`, `n_selected_features`, `clusters[]` (sizes), `embeddings`,
`parameters`, `figures`.

## peaks

`call_peaks` runs `setup → blacklist → load → qc-metrics → filter → features →
embed → cluster → peaks → save` (no qc-plots). `peaks.py`:
`snap.tl.macs3(adata, groupby=group, qvalue=peak_qvalue, call_broad_peaks,
blacklist, inplace=True, n_jobs)` stores `adata.uns["macs3"]` as
`dict[str, polars.DataFrame]` (returns None); `snap.tl.merge_peaks(tables,
chrom_sizes, half_width=peak_half_width)` returns a polars DataFrame whose
`Peaks` column is `"chr:start-end"` strings; the merged BED is written and
`snap.pp.make_peak_matrix(adata, use_rep=merged["Peaks"].to_list(),
inplace=False, file=tables_dir/peak_matrix.h5ad)` builds the separate backed
peak object (obs copied). `use_rep` must be a STRING uns-key or a LIST of peak
strings — a DataFrame PANICS. The whole block is try/except: a peak-matrix input
without fragments reports a gap and continues. Reports `peaks_per_cluster`,
`merged_peaks`, the calling parameters.

## markers

`find_markers` runs the peak chain + `markers`. `markers.py` targets
`peak_adata if not None else adata`. `snap.tl.marker_regions(target, groupby,
pvalue=marker_pvalue)` → `dict[str, list[str]]` (cluster-vs-rest, FAST), written
as one `markers_<cluster>.txt` per cluster. Pairwise `snap.tl.diff_test(target,
cell_group1=<indices>, cell_group2=<indices>, min_log_fc, min_pct,
direction="both")` is **opt-in via `run_diff`** (default false): it is
O(clusters²) logistic regressions (~2 min/pair on a real peak matrix) and takes
cell INDICES or barcode NAMES, not cluster labels — the fragment selects indices
from the group column with `np.flatnonzero`. `diff_results` maps
`"<c1>_vs_<c2>"` → region count, CSVs written when non-empty. Reports
`total_markers`, `groups`, `tables`, `diff_accessibility`, and the descriptive
(double-dipping) caveat.

## annotate

`annotate_clusters` runs `setup → blacklist → load → qc-metrics → filter →
features → embed → cluster → gene-activity → annotate → save` (no peaks/markers).
- `gene-activity.py`: `snap.pp.make_gene_matrix(adata, gene_anno, inplace=False,
  file=tables_dir/gene_activity.h5ad, id_type="gene", upstream=gene_upstream,
  include_gene_body=gene_body, gene_name_key="gene_name")` → separate backed
  `gene_adata`. Skipped with a gap note without `gene_anno` (or without
  fragments).
- `annotate.py`: three routes in priority order — (1) `mapping_path` explicit
  cluster→label JSON wins; (2) panel scoring on `gene_adata.to_memory()`
  (normalize_total + log1p, per-cluster panel mean minus rest mean, `min_score`
  gates acceptance, generic example panel when `panel_path` null); (3)
  `reference_path` scRNA label transfer (shared genes ≥50, concat + pca +
  sklearn KNeighborsClassifier cosine, per-cluster mode vote). Writes
  `obs[label_column]`; a script never invents labels — without gene activity,
  reference or mapping it reports the gap. Reports `method`, `assignments`,
  `scores`, `gene_activity_present`.

## motif

`motif_enrichment` runs the marker chain + `motif`. `motif.py` needs `markers`
(from the markers step), `genome_fasta` (local; the snapatac2 lazy download is
~845 MB) and a motif DB (`snap.read_motifs(motif_db)` MEME → list[PyDNAMotif],
else `snap.datasets.cis_bp(unique=True)`, a download).
`snap.tl.motif_enrichment(motifs, regions, genome_fasta, method=motif_method)` →
`dict[str, polars.DataFrame]`, one `motif_<cluster>.csv` per cluster. This is
motif enrichment (sequence prediction), NOT chromVAR per-cell deviations
(chromVAR is R/pip, a project-env gap). Missing any input → gap note, continue.
Reports `results` (per-cluster enriched-motif counts), `motif_method`,
`motif_pvalue`.

## integrate

`integrate_samples` runs the clustered chain + `integrate`. `integrate.py` needs
a `batch_column` with ≥2 batches and a computed `X_spectral`.
`snap.pp.harmony` / `scanorama_integrate` / `mnc_correct` (batch=batch_column,
use_rep="X_spectral") → `obsm["X_spectral_harmony"]` / `_scanorama` / `_mnn`.
Corrected values are for joint clustering/viz only — never peak calling or diff.
Reports `integration` (method/batches/embedding), `embeddings`.

## env_probe

Probes a 23-package catalogue across 6 modules via `importlib.util.find_spec` +
`importlib.metadata.version`, reports present/absent as facts, and `write_panel`
emits a starter marker-panel JSON (the same generic immune/stromal panel as
`annotate.py`). Core required: snapatac2/scanpy/anndata/numpy/pandas/polars/
pyarrow/scipy/h5py/scikit-learn/leidenalg/igraph/umap-learn/macs3/plotly/
matplotlib; kaleido core optional. Gaps (modules=[]): episcanpy (numpy 2
incompatible), chromVAR/ArchR/Signac (R), scvi-tools/magic-impute (pip-only).
The one declared tool output schema is `schemas/env-probe.json` (draft 2020-12,
additionalProperties:false). Absent packages are project-environment next steps
(`env_request`), never in-place installs (the plugin env is read-only,
conda-locked, `refusePipDependencies`).

## snapatac2 API mapping (quick reference)

| step | snapatac2 call | result |
| --- | --- | --- |
| load | `pp.import_fragments` / `pp.make_fragment_file` / `read` / `read_10x_mtx` | backed adata, obs n_fragment/frac_dup/frac_mito |
| qc-metrics | `metrics.frag_size_distr`, `metrics.tsse` | obs tsse; uns frag-size |
| filter | `pp.filter_cells(min_counts,min_tsse,max_counts,max_tsse)` | in-place subset |
| features | `pp.add_tile_matrix`, `pp.select_features`, `pp.scrublet`, `pp.filter_doublets` | var selected; obs doublet_probability |
| embed | `tl.spectral`, `pp.knn` | obsm X_spectral; obsp connectivities |
| cluster | `tl.leiden`, `tl.umap`, `pl.umap` | obs leiden[_r]; obsm X_umap |
| peaks | `tl.macs3`, `tl.merge_peaks`, `pp.make_peak_matrix` | uns macs3; separate peak_adata |
| markers | `tl.marker_regions`, `tl.diff_test` | dict cluster→regions; polars diff DataFrame |
| gene-activity | `pp.make_gene_matrix` | separate gene_adata |
| motif | `read_motifs`, `datasets.cis_bp`, `tl.motif_enrichment` | dict cluster→polars DataFrame |
| integrate | `pp.harmony`, `pp.scanorama_integrate`, `pp.mnc_correct` | obsm X_spectral_<method> |
| save | `adata.write`, `.close()` | exported .h5ad |

`transfer_labels` is WITHIN-object kNN propagation, NOT cross-modal ATAC↔RNA
transfer — annotation uses gene-activity space + sklearn kNN instead. There is
NO chromVAR in snapatac2; the TF layer is `motif_enrichment` only.

## Preset mapping

| Module | Preset | Chain coverage |
| --- | --- | --- |
| 01 scATAC QC | `module: 1` | setup → load → qc-metrics → qc-plots → filter → save |
| 02 Preprocess & cluster | `module: 2` | + blacklist/features/embed/cluster |
| 03 Peaks & marker regions | `module: 3` | + peaks/markers |
| 04 Gene activity & annotation | `module: 4` | + gene-activity/annotate |
| 05 Motif & integration | `module: 5` | + motif/integrate |
| 06 Full scATAC chain | `module: 6` | every capability |

A `module` preset uses its explicit step list with NO dependency closure, so
blacklist is listed explicitly where needed; `steps` closes over `requires`.
Steps needing a missing genome resource, reference or package report the gap and
continue (the notebook footer lists what did not run). Raw reads→fragments
(Cell Ranger ATAC, chromap) are wrappers/upstream, not notebook cells. Bulk
ChIP/CUT&Tag is the sibling `chromatin-profiling` plugin (R, pure downstream).

## JSON contract

Every tool emits one object with `input`/`output` plus its section blocks.
Common conventions: `figures[]` lists HTML paths actually written (snapatac2
plots are plotly Figures, `show=False`, written as interactive `.html`; a
kaleido PNG is optional via `static_png`); `notes[]` carries warnings and
skipped-step reasons; a skipped step appends a "… skipped: …" note, never an
empty success. Only `env_probe` declares a tool-level output schema
(`schemas/env-probe.json`); the other tools' contracts are documented here and
in the playbook. Adapter and notebook share the fragment library, so the same
reading fed to both produces the same numbers by construction. The renderer
(`build_notebook.py`, RENDERER_VERSION=3) reads each fragment plus
`references/step-notes/<step>.md` (a missing note → empty Markdown cell, code
still runs), builds one config cell from `steps.json` parameters + flags/config
JSON (precedence flag > config > default), and ends with a close cell.