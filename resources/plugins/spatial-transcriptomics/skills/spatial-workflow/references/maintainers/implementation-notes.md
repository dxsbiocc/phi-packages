# Implementation notes — spatial Python route

What each tool computes and how the layers divide the work. The chain is
deterministic: every step reads one object, writes one, and reports the
parameters it used. The analysis lives in `assets/fragments/` (per-step teaching
notes in `references/step-notes/`) — the tools are thin adapters over the same
fragments the notebook renders, so **adapt the fragment, not the script**.

## Architecture

```
assets/steps.json            registry: parameters (name/type/default/comment),
                             steps (id/title/requires), presets, kernel, footer
assets/fragments/<step>.py   self-contained teaching code, exec'd in a shared
                             namespace; the SINGLE analysis source for both paths
references/step-notes/<step>.md  the notebook's per-step Markdown teaching cell
scripts/_common.py           flag parsing, read/write object, emit/fail/info,
                             detect_kind, sorted_levels, figure helpers
scripts/_config.py           config_namespace(): builds EVERY parameter global
                             from steps.json defaults + flag overrides; structural()
                             adds data_path/output/emit/tool_input/tool_output
scripts/_fragments.py        run_fragments(steps, ns) execs fragments with stdout
                             redirected to stderr; run_report(name, ns) execs reports/
scripts/readers/             format readers (visium, visium_hd, stereoseq) returning
                             (adata, meta); sniff_h5 classifies HDF5 layouts
scripts/<tool>.py            thin adapter: parse flags → build namespace →
                             run_fragments → run_report; stdout = exactly one JSON
scripts/reports/<name>.py    serializes fragment bookkeeping into the JSON contract
build_notebook.py            renderer (RENDERER_VERSION=3): reads fragments +
                             step-notes, emits config cell + code/Markdown cells
```

**Why `_config.config_namespace`.** A fragment references ~50 parameter globals
(`min_genes`, `neighbor_method`, `svg_method`, …). Building the namespace from
the registry guarantees every global exists with its declared default, coerced
to the declared type, then overridden by flags (flag name == parameter name) and
structural values. This is what keeps the adapter namespace a superset of the
notebook config namespace — the two paths agree by construction, and a fragment
never hits a `KeyError` because a flag was omitted.

**stdout discipline.** `run_fragments` redirects fragment stdout (print, scanpy
logs) to stderr; the report's `emit` writes exactly one JSON object to stdout.
This mirrors the R adapters' `sink(stderr())`. A tool's stdout is always
parseable JSON; progress is on stderr.

**Bookkeeping.** Fragments leave `figure_paths`, `analysis_steps`,
`analysis_notes`, `qc_filters`, `qc_scored`, `obj_before`, `pca_variance`,
`cluster_columns`, `embeddings`, `markers`, `marker_tables`, `annotation_*`,
`svg_results`, `nhood_results`, `cooc_results`, `spatial_diff_results`,
`deconv_results`, `ccc_results`, `platform_info`, `bin_info` in the namespace.
Reports read these; they never recompute analysis.

## Readers

`detect_kind` + `sniff_h5` classify the input; `read_spatial` dispatches:

| kind | reader | notes |
| --- | --- | --- |
| visium-dir / visium-h5 / 10x-h5 / 10x-matrix-dir | `visium` | Space Ranger matrix + `spatial/` (tissue_positions, scalefactors, hires image) |
| stereoseq-gef / stereoseq-cellbin | `stereoseq.read_gef` | HDF5 `geneExp` (square bin) or `cellBin` (segmented cells); `bin_size` selects/reads the bin |
| stereoseq-gem | `stereoseq.read_gem` | TSV gene-exp matrix; DNB coords aggregated by `bin_coordinates` |
| dense-table | `visium_hd.read_table` | Visium HD square-bin/seg .csv/.tsv |
| anndata / loom | `visium` | `.h5ad`/`.loom` with `obsm["spatial"]` |

`finalize` makes var names unique, asserts `obsm["spatial"]` exists (spatial
statistics need coordinates), and records `n_units`/`n_genes`/`coord_key` in
meta. The caller adds the `counts` layer. `bin_coordinates` aggregates DNB
centres to square-bin origins (`(x // N) * N`). Stereo-seq cellbin segmentation
is produced upstream by SAW cellCut; the reader consumes it, never re-segments.
Where a file's internal layout is undocumented the reader raises rather than
guessing field names.

## py_inspect

Dispatches on input kind, snapshots `X` as `counts` when absent, reports
`platform`/`bin` provenance, `shape` (units/genes/layers/obsm/obsp/uns/snapshot),
`counts_present`, `coordinates_present`, `value_class`, `metadata` columns worth
using, and `candidate_steps` computed from what the object already contains
(counts → qc/preprocess/deconvolve; coords → spatial; leiden → markers/annotate;
lognorm+leiden → spatial-statistics/communication). `min_genes`/`min_units`
default 0 here — gating belongs to py_qc. A `bin_auto` note flags an
auto-selected Stereo-seq bin so the agent overrides it deliberately.

## py_qc

`calculate_qc_metrics` on counts with mitochondrial/ribosomal/hemoglobin regex
sets; the threshold filters, each reported with its removed-unit count.
`distribution` (min/q01/median/q99/max) is the evidence for choosing thresholds.
Spatial QC adds tissue-coverage maps (`qc-spatial-before/after`): low-quality
units often sit at tissue edges. The mt gate is guarded — if no gene matched
`mt_pattern` the column is absent and the gate is skipped with a warning (mouse
`^mt-` vs human `^MT-`), never a `KeyError`. `kept_fraction` < 0.5 warns.

## py_preprocess

`normalize_total` + `log1p`; `highly_variable_genes` `flavor=seurat_v3` on the
**counts** layer (needs scikit-misc, shipped); `layers["lognorm"]` snapshot
**before** scaling; `scale(max_value=10)`; `pca(mask_var="highly_variable")` —
the object is never gene-subset. Leiden at each requested resolution, UMAP.
Then the **spatial adjacency graph** (`spatial-neighbors`): mode-specific
squidpy builder (`spatial_neighbors_knn/radius/delaunay/grid`; the catch-all
`spatial_neighbors` is deprecated in 1.8) on `obsm[coord_key]`, stored in
`obsp["spatial_connectivities"]`. `spatial-plots` maps clusters/markers in tissue
space via `sc.pl.embedding` (works without a tissue image). Reports
`spatial_graph` (method/edges/mean_degree), `variance`, `clusters[]`.

## py_markers

`rank_genes_groups` (`test_use=wilcoxon` default) on `layers.lognorm` when
present. Tables filtered by `logfc_threshold`/`min_pct`/`only_positive`, sorted
by adjusted p-value, one CSV per group. `total_markers=0` early-exits with a
note. The report flags whether the test ran on the lognorm snapshot and that
p-values are descriptive (double dipping).

## py_annotate

Per-cluster panel mean minus the mean of all other genes, on the pre-scale
snapshot — difference in means keeps panel sizes comparable. `min_score` gates
acceptance; below it a cluster stays unassigned. `mapping_path` overrides
scoring but must label every cluster. Runs `annotate` only with `panel_path`,
and `assign-labels` with `panel_path` or `mapping_path`; with neither, nothing
is annotated (reported). Writes `label_column` and maps labels in space.

## py_spatial

Runs `spatial-neighbors` then `svg`, `nhood`, `spatial-diff`:

- **svg**: `sq.gr.spatial_autocorr` (mode moran/geary) on lognorm, analytic
  p-values (`n_perms=None`), corrected with a statsmodels method (`fdr_bh`;
  `benjamini-hochberg`/`bh` aliased to `fdr_bh`). Reads `uns["moranI"]`/
  `["gearyC"]`, already sorted; reports `top[]` + `n_significant`.
- **nhood**: `sq.gr.nhood_enrichment` (permutation z-score per cluster pair) +
  `sq.gr.co_occurrence` (distance-resolved). Needs ≥2 groups; reports
  `top_enriched_pairs[]` and co-occurrence `distance_range`/`shape`.
- **spatial-diff**: with `group_column`, subsets each condition, rebuilds its
  OWN knn graph, recomputes Moran's I, ranks genes by cross-condition range,
  reports condition-specific SVGs and cluster `composition`. Needs ≥2 usable
  conditions (`min_units`); a descriptive screen, not a significance test.

Each skipped section carries `{"skipped": reason}` so absence is never read as
a clean result.

## py_deconvolve

scipy NNLS, dependency-free (no torch). Builds a pseudobulk reference signature
(mean library-size-normalized profile per cell type) over up to 500 shared
signature genes (spatial HVGs preferred, ranked by spatial variance); for each
unit solves non-negative weights minimising `||S w - y||`, rescaled to sum 1.
Writes `obsm["deconvolution"]` and `obs["<type>_prop"]`. Without `reference_path`
it skips with a message; `deconv_method` must be `nnls` (cell2location/tangram
are project-env gaps). Prefers the reference `counts` layer; without it assumes
log1p X and inverts with `expm1` (noted). No shared genes → identifier-mismatch
error. Proportions are estimates, not cell counts.

## py_communication

`squidpy.gr.ligrec` (CellPhoneDB-style permutation test). Offline by default: a
curated inline ligand-receptor DataFrame (~70 pairs, source/target) ships in the
fragment; `interactions_path` (TSV/CSV with source/target or ligand/receptor)
overrides it. ligrec reads `adata.X` with `use_raw=False`, so the adapter builds
a minimal AnnData whose X is the lognorm snapshot (no layer argument exists).
Needs a categorical `cluster_key` and ≥2 groups. Reports `n_significant`,
`top[]` (ligand/receptor/source_cluster/target_cluster/mean). Symbols match
case-insensitively (ligrec uppercases). omnipath is conda-locked but the online
fetch is avoided; liana/CellPhoneDB are project-env gaps.

## py_env_probe

Probes a 22-package catalogue across 6 modules (core, spatial, deconvolution,
communication, …), reports present/absent as facts, and `write_panel` emits a
starter marker-panel JSON. Absent packages are project-environment next steps,
never in-place installs (the plugin env is read-only and conda-locked).

## Preset mapping

| Module | Preset | Chain coverage |
| --- | --- | --- |
| 01 Spatial QC | `module: 1` | setup → load → qc-metrics → qc-plots → filter → save |
| 02 Preprocess & cluster | `module: 2` | + normalize/hvg/scale/pca/cluster |
| 03 Spatial domains & markers | `module: 3` | + spatial-neighbors/spatial-plots/markers/annotate/assign-labels |
| 04 SVG & neighbourhoods | `module: 4` | + svg/nhood |
| 05 Deconvolution & communication | `module: 5` | + deconvolve/communication |
| 06 Full spatial chain | `module: 6` | + spatial-diff (every capability) |

Steps needing a missing reference or package report the gap and continue (the
notebook footer lists what did not run). Raw reads→counts (Space Ranger, SAW)
and Stereo-seq cell segmentation are wrappers/upstream, not notebook cells.

## JSON contract

Every tool emits one object with `input`/`output` plus its section blocks.
Common conventions: `figures[]` lists PNG paths actually written; `notes[]`
carries warnings and skipped-step reasons; a skipped advanced section is
`{"skipped": reason}`, never an empty success. The notebook's `phi_render`
metadata (parameter sources, fragment sha256, data fingerprint) is the parallel
audit record. Adapter and notebook share the fragment library, so the same
reading fed to both produces the same numbers by construction.