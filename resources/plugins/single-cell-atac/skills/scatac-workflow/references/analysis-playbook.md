# Single-cell ATAC analysis playbook (Python route)

The agent's operating manual. Each stage below is a decision procedure:
**Call** the tool, **Read** the named JSON fields and figures, **Decide** the
parameters from what you read, and respect the **Red lines**. The one-line
*why* gives the principle; the hard rules (fragments first, renderer/agent
split, one reading for both paths, tiles-vs-peaks, backed-object close) live in
SKILL.md `Rules`; the learner-facing short version lives in
`references/step-notes/` (one note per step). Sources and licenses are at the end.

Quick reference: [`red-line-cheatsheet.md`](red-line-cheatsheet.md) tabulates the
hard red lines and common misreadings for §1–§10; this playbook remains the full
Call/Read/Decide procedure.

Every number here is a **starting point for scATAC data, not a default to
accept**. Read the distributions first, adapt, and say in your answer which
value you chose and why. Quote the tool's JSON numbers instead of re-deriving
them. The engine is snapatac2 2.10 (file-backed, Rust); scanpy 1.12 is used only
for gene-activity downstream and matplotlib panels.

## 0. The decision loop

Chain: `scatac_inspect` → `scatac_qc` → `scatac_preprocess` → then branch on the
question:

- **What are the cell states / which regions define each cluster?** →
  `scatac_peaks` → `scatac_markers`.
- **What cell type is each cluster?** → `scatac_annotate` (gene activity + panel
  or scRNA reference).
- **Which TFs may drive each cluster's open chromatin?** → `scatac_motif`
  (needs a genome FASTA).
- **Do samples cluster by batch instead of biology?** → `scatac_integrate`
  (needs a batch column with ≥2 samples).
- **Which regions differ between two specific clusters / conditions?** →
  `scatac_markers` with `run_diff=true` (pairwise `diff_test`, opt-in, slow).

`scatac_notebook` is the parallel teaching path; feed it the *same* threshold
readings (flags or `config` JSON). Each tool re-runs its chain from the raw input
you give it (snapatac2 objects are file-backed and heavy to hand off mid-chain),
so pass the same input and the same thresholds every time to keep the two paths
and successive calls consistent. Skip a step only when `scatac_inspect`'s
`candidate_steps` says its result already exists.

Loop at every stage: call with loose/starting values → read the JSON → decide →
re-call with the decided values → verify `after`, `warnings`, `notes`. A default
that reached the output untouched is a decision you did not make.

Call `scatac_env_probe` before promising motif enrichment, gene activity,
integration or any optional package: the plugin environment is read-only; a gap
is reported as a fact with a project-environment next step (`env_request`),
never an in-place install. chromVAR / ArchR / Signac (R) and episcanpy (numpy 2
incompatible) are gaps, not capabilities — say so, do not impersonate them.

## 1. Inspect the input — `scatac_inspect`

**Call.** `scatac_inspect(input, output=<project>/loaded.h5ad)`. Keep
`min_num_fragments` at its import default here: per-cell gating belongs to
`scatac_qc`, which reports per-gate removals; dropping barcodes silently at load
hides the empty-droplet shoulder. Choose `genome` (hg38/mm10/hg19/mm39 bundled
offline) or pass `chrom_sizes` for another build — the genome decides every
later blacklist, TSSe and motif lookup. Pass `gene_anno` only if you want TSSe
in the same call; otherwise it is skipped with a gap note.

**Read.** `kind` (fragment / fragment-like-table / bam / anndata /
10x-matrix-dir / 10x-h5 / peak-file); `genome` (`source`: bundled / file /
saved-object / 10x-matrix, `genome`, `n_chrom`, `input_kind`); `shape`
(`cells`, `features` — fragments import with **0 features**, the tile matrix
comes later); `metadata`/`metadata_columns` (which of `n_fragment`, `frac_dup`,
`frac_mito`, `tsse`, `sample`, `condition`, `leiden`, `cell_type` are present);
`candidate_steps`; `notes`.

**Decide.**
- **Input kind and genome first.** A `fragment`/`bam` input has raw accessibility
  and can run the whole chain. A `10x-matrix-dir`/`10x-h5`/`peak-file` input
  already carries features but **no fragments** — peak calling (MACS3) and
  scrublet need fragments, so those steps will report a gap; say which input you
  analysed and what it cannot support.
- `sample_column`/`batch_column` for integration and per-sample QC: the
  donor/condition column, not a cluster column.
- Whether `gene_anno` is available (TSSe, gene activity) and `genome_fasta`
  (motif); if not, name the gap and the project-directory download / `input-path`
  next step rather than triggering a slow lazy download into the read-only env.

**Red lines.** No peak calling or scrublet without fragments. No TSSe or gene
activity without `gene_anno`. No motif without `genome_fasta`.
*Why:* MACS3 and scrublet read the fragment/feature evidence; TSSe and gene
activity are defined on gene coordinates; motif scanning needs the sequence.

## 2. Quality control — `scatac_qc`

**Call 1 (reading).** Leave `min_tsse`/`max_tsse` unset (a None bound skips that
gate) and pass a loose `min_counts` so you read the distribution before gating
it. Pass `gene_anno` to enable TSSe — the best ATAC signal-to-noise metric. Set
`figures_dir` and look at the before-panels: the fragment-size ladder
(`qc-frag-size-before`), TSSe vs fragment count (`qc-tsse-before`) and the
log-x metric histograms (`qc-metrics-before`).

**Read.** `distribution[]` for `n_fragment`, `frac_dup`, `frac_mito`, `tsse`:
`min/q01/median/q99/max`. `features_scored` (which metrics actually computed —
`tsse` absent means no `gene_anno`). `before.per_sample` when a sample column
exists. `warnings` (a "no tsse computed" warning is a missing annotation, not a
clean result).

**Decide.**
- `min_counts`: from the `n_fragment` distribution — read the valley between the
  empty-droplet shoulder and the cell peak on the log-x histogram. scATAC cells
  are shallow; start ~1000 and adapt, do not accept the default blindly.
- `max_counts`: the upper shoulder can be doublets/aggregates; set it only when
  the distribution shows a distinct high tail, otherwise leave None and let
  scrublet (in `scatac_preprocess`) handle doublets.
- `min_tsse`: ATAC-specific signal-to-noise. A cell with enough fragments but low
  TSSe is still noise. Read the `tsse` shoulder before fixing the gate; if `tsse`
  is absent, fix `gene_anno` first — a gate that never fires is not QC.
- `frac_dup`/`frac_mito` are reported for context (over-sequencing, damage); they
  are not universal gates the way `n_fragment` and `tsse` are.

**Call 2 (gating).** Pass the decided values. **Read** `filters[]` (`filter`,
`removed` — trace a surprising total to one gate), `after.kept_fraction`,
`after.removed_cells`, `warnings`, and the after-panels. Doublet removal is NOT
here: scrublet needs the selected feature matrix and runs in `scatac_preprocess`,
so this module stays lightweight (no tile matrix).

**Red lines.** `kept_fraction` < 0.5 → stop and re-read the distributions; zero
cells kept → relax or check the genome build matches the fragments. QC runs on
raw fragments/counts only. Samples that differ in depth get per-sample gates.
*Why:* fragment count, duplication, mitochondrial fraction and TSSe carry
distinct signals; TSSe is the one that separates real ATAC signal from noise.

## 3. Tile matrix, embedding, clustering — `scatac_preprocess`

**Call.** Pass the thresholds decided in §2 (the same numbers), plus `bin_size`,
`n_features`, `n_comps`, `n_neighbors`, `resolution` (a list to sweep),
`doublet_filter`/`doublet_threshold`/`expected_doublet_rate`, and `blacklist`
(null uses the shipped ENCODE list for hg38/mm10). Set `figures_dir` for the
UMAP cluster plot.

**Read.** `steps` (tile-matrix / select-features / doublet-filter / spectral /
knn / clusters — which actually ran); `n_selected_features`; `clusters[]`
(`column`, `clusters`, `sizes[{cluster,cells}]`); `embeddings` (`X_spectral`,
`X_umap`); `parameters`; `figures` (`dimplot-clusters`); `notes`. The
`filters[]`-style doublet count is in `notes`/stderr; the scrublet removal is
also recorded as a gate.

**Decide.**
- **Tiles for clustering, peaks for comparison.** `bin_size` is a resolution
  decision: ~500 bp tiles cover the genome uniformly for the embedding; called
  peaks (§4) are the features of record for differential accessibility. Say what
  bin you chose and why. A 10x matrix input already has features, so tiling is
  skipped.
- `n_features`: how many tiles enter the spectral embedding (default 500000);
  `select_features` drops uninformative bins and blacklist artifacts. Lower it
  for speed on small data, raise it to retain rare regulatory signal.
- `n_comps`: spectral embedding has **no variance-ratio scree plot** (unlike
  PCA) — choose from downstream clustering stability, not a curve. Start 30.
- `resolution`: sweep it (e.g. `[0.5, 1.0, 1.5]`) and compare cluster counts and
  sizes rather than trusting one value. One resolution writes `obs['leiden']`;
  several write `obs['leiden_<r>']`.
- `doublet_filter`: keep true unless you have a reason; `expected_doublet_rate`
  scales with loaded cell number. A doublet inflates accessibility and fakes a
  "high-activity" new type.

**Red lines.** The kNN graph and Leiden run on `X_spectral`, never on UMAP (UMAP
distorts distances; it is for looking). Corrected/integrated values never enter
the embedding used for peak calling or diff. Blacklist artifacts must be removed
before the embedding or they dominate it with technical signal.
*Why:* spectral embedding denoises the sparse binary accessibility matrix into a
space where neighbours are meaningful; UMAP and batch-corrected values are not
that space.

## 4. Per-cluster peaks — `scatac_peaks`

**Call.** Pass the clustered chain plus `group_by` (default leiden),
`peak_qvalue` (default 0.05), `peak_half_width` (default 250 bp),
`call_broad_peaks` (false for ATAC narrow peaks; true only for histone marks),
and `tables_dir` for the merged BED and peak-count matrix. MACS3 needs fragments.

**Read.** `peaks_per_cluster` (`{cluster: n_peaks}`); `merged_peaks` (the
non-overlapping merged set size); `group_by`, `peak_qvalue`, `peak_half_width`;
`steps`; `notes`. The merged BED (`merged_peaks.bed`) and the separate backed
peak matrix (`peak_matrix.h5ad`) are written to `tables_dir`.

**Decide.**
- **Per-cluster calling is a biological decision.** Calling peaks on each
  cluster's fragments captures cell-state-specific open chromatin that a global
  call misses. `peak_half_width` extends each summit before merging — wider
  merges nearby regulatory elements, narrower keeps them distinct; say what you
  chose.
- `peak_qvalue`: relax it if `merged_peaks` is 0 or a cluster yields very few
  peaks (sparse clusters need a permissive cutoff); tighten it to keep only
  strong summits. A `merged_peaks: 0` with a "no merged peak set" note means
  peak calling did not produce features — check fragments, cluster column and
  qvalue, do not read it as "no open chromatin".
- The peak matrix is a **separate object** (`peak_adata`); the backed tile matrix
  cannot be resized in place. Tiles stay the embedding, peaks drive markers.

**Red lines.** Peak calling runs on raw fragments, never on corrected/scaled
values. A peak-matrix input without fragments reports a gap here — do not fake
peaks. ATAC uses narrow peaks; broad peaks are for histone marks.
*Why:* MACS3 models the fragment insertion distribution; called peaks are the
open-chromatin features of record for comparison.

## 5. Marker regions & differential accessibility — `scatac_markers`

**Call.** Pass the peak chain plus `group_by`, `marker_pvalue` (default 0.01),
`tables_dir`, and `run_diff` (default **false**). Only set `run_diff=true` when a
specific cluster-vs-cluster contrast is needed — pairwise `diff_test` is
O(clusters²) logistic regressions and is slow. `diff_min_log_fc` (0.25) and
`diff_min_pct` (0.05) gate the diff.

**Read.** `total_markers`; `groups[]` (`{cluster, n_markers}`); `tables[]`
(`{cluster, n_markers, path}` — one `markers_<cluster>.txt` per cluster);
`diff_accessibility` (`{ "<c1>_vs_<c2>": n_regions }`, empty unless `run_diff`);
`marker_pvalue`, `diff_min_log_fc`, `diff_min_pct`; `notes`.

**Decide.**
- `marker_regions` is always **cluster-vs-rest** and fast; it gives each
  cluster's specifically accessible peaks. Read the per-cluster tables and map
  the regions to genes/anatomy to form a hypothesis.
- `marker_pvalue`: relax it if `total_markers` is 0 (clusters may not differ, or
  the peak matrix is sparse); the p-values are **descriptive** — clusters were
  defined from the same data (double dipping), so treat the ranking as evidence,
  not significance.
- `run_diff`: enable only for a named contrast. `diff_test` takes cell indices
  (the fragment selects them from the group column), returns log fold change and
  significance per region; it is the scATAC analogue of differential expression.
  Both marker and diff run on the peak matrix when peaks were called, else tiles.

**Red lines.** Differential accessibility runs on raw peak/tile counts, never on
integrated values. A `total_markers: 0` or empty `diff_accessibility` with a note
is a gap or a null result, not "no regulation". Do not claim significance from
double-dipped cluster p-values.
*Why:* the test assumes the original count distribution; corrected values break
the mean-variance relationship, and same-data clusters invalidate naive p-values.

## 6. Gene activity & annotation — `scatac_annotate`

**Call.** Pass the clustered chain plus `gene_anno` (for the gene activity
matrix), `gene_upstream` (2000 bp promoter), `gene_body` (true), and ONE
annotation route: `mapping_path` (explicit cluster→label, wins), `panel_path`
(cell type→gene symbols; null uses a generic example panel), or `reference_path`
+ `reference_label_key` (scRNA label transfer). `label_column` (default
cell_type), `min_score`, `tables_dir`, `figures_dir`.

**Read.** `method` (mapping / panel-score / reference-transfer / null);
`gene_activity_present`; `panel[]`; `assignments[]` (`{cluster, label, cells}`);
`scores[]` (`{cluster, cells, best, scores[{cell_type, score, panel_mean,
genes_found, genes_missing}]}`); `label_column`; `notes`. The gene activity
matrix is a separate `gene_activity.h5ad` in `tables_dir`.

**Decide.**
- **Gene activity is a proxy, not RNA.** It aggregates promoter/gene-body
  accessibility into a gene-by-cell matrix so scRNA marker panels and label
  transfer work in a shared gene space. Use it for annotation and marker
  visualization, never for differential-expression claims.
- Route priority: an explicit `mapping_path` is a curator decision and wins;
  panel scoring is a heuristic ranking (read `scores`, `genes_found` — a panel
  gene missing from gene activity means the wrong identifier type or a sparse
  locus); reference transfer needs a scRNA reference from the single-cell plugin
  and ≥50 shared genes.
- `min_score`: below it a cluster stays `unassigned` — do not force a label.
  Edit the panel to tissue-specific markers for real data; the generic panel is
  an example, not a universal immune/stromal truth.

**Red lines.** No annotation without `gene_anno` (gene activity), a
`reference_path`, or a `mapping_path` — the step reports the gap and a script
must not invent labels. No label-transfer claim without a scRNA reference;
without one, only gene-activity marker inference, labelled as such.
*Why:* annotation is inference from accessibility, not observation of RNA; the
label source must exist and be stated.

## 7. Motif enrichment — `scatac_motif`

**Call.** Pass the marker chain plus `genome_fasta` (a local FASTA — the
snapatac2 lazy download is ~845 MB from a slow host), `motif_db` (.meme; null
uses snapatac2 `cis_bp`, itself a download), `motif_method` (binomial /
hypergeometric), `motif_pvalue` (1e-5), `group_by`, `marker_pvalue`, `tables_dir`.

**Read.** `results` (`{cluster: n_enriched_motifs}`); `motif_method`,
`motif_pvalue`; `steps`; `notes`. Per-cluster enrichment CSVs (`motif_<cluster>.csv`)
are written to `tables_dir`.

**Decide.**
- Motif enrichment asks which TF motifs are over-represented in each cluster's
  marker regions versus background — the regulatory hypothesis layer. Report
  which database (JASPAR/cis_bp) and background you used.
- `motif_method`: binomial (default) vs hypergeometric; `motif_pvalue` controls
  stringency. Read the top motifs per cluster and cross-check with the marker
  genes / annotation.
- Without `genome_fasta` the step reports the gap and continues — name the
  project-directory FASTA next step, do not download into the read-only env.

**Red lines.** Motif hits are **sequence-level predictions**, not measured TF
binding; they require experimental validation. This is motif enrichment, NOT
chromVAR per-cell TF deviations — chromVAR is an R/pip gap, not shipped; do not
impersonate it. A full TF-gene GRN is a snapatac2 extension via `inject`.
*Why:* enrichment tests motif frequency in accessible sequence; it does not
measure occupancy or per-cell activity.

## 8. Sample integration — `scatac_integrate`

**Call.** Pass the clustered chain plus `batch_column` (default sample),
`integrate_method` (harmony / scanorama / mnc), `n_comps`, `n_neighbors`,
`resolution`. Needs ≥2 batches and a computed `X_spectral`.

**Read.** `integration` (`{method, batch_column, n_batches, batches[],
embedding}` — `X_spectral_harmony` / `_scanorama` / `_mnn`); `embeddings`;
`steps`; `notes`.

**Decide.**
- Integration removes technical variation so clusters reflect biology, not run.
  Harmony (default) corrects the spectral embedding; scanorama/mnc are
  alternatives. Re-cluster on the corrected embedding if batch drove clusters.
- **The sample is the experimental unit.** Fewer than two biological replicates
  per condition means no condition-level claim; cells are pseudoreplicates.
- Do not over-correct: integration can erase real biology. Check that known
  cell-type markers still separate after correction.

**Red lines.** Corrected/integrated values are for joint clustering and
visualization ONLY — they never enter peak calling or differential
accessibility, which return to raw counts and model by sample. No integration
without a batch column with ≥2 samples (the step reports the gap).
*Why:* batch correction changes the mean-variance structure the tests assume;
condition-level inference needs biological replicates, not cells.

## 9. The notebook path — `scatac_notebook`

**Call.** `scatac_notebook(output=<project>/scatac.ipynb, data=<input>,
module=<1-6> | steps=[...], ...)`. Pass the *same* thresholds read in §2 as flags
or a `config` JSON (keys mirror `assets/steps.json` parameters). Pick a preset:
1 QC, 2 preprocess+cluster, 3 peaks+markers, 4 gene activity+annotation, 5
motif+integration, 6 full chain; or name `steps` (dependencies are added
automatically — but a `module` preset uses its explicit list, so blacklist is
already included where needed).

**Read.** The rendered notebook: one config cell at the top holding every
tunable value with a source tag, then analysis cells referencing those
variables, then a close cell. `self_check=true` executes every cell in the
output directory before delivery and fails loudly if one breaks.

**Decide.**
- `inject` adds an external Python script verbatim as a code cell (e.g. a TF-gene
  GRN extension); `merge` splices in another notebook whole or a cell range
  (`path.ipynb[0:6]`). Sources must be Python — one kernel per notebook.
- The renderer fills placeholders from the data and the requested steps; it does
  NOT choose steps or thresholds or interpret the request — those are yours.
- Render into the user's working directory, never copy the template dir. The
  notebook ends by closing the backed object.

**Red lines.** One reading, both paths: the script chain and the notebook must
agree by construction (same flags/config), not by luck. No install cells, no
Colab paths, no raw-read alignment cells (that is an upstream wrapper). A step
that needed a missing resource says so and continues — its absence is not a
clean result.
*Why:* the notebook is the teaching twin of the script chain; divergence between
them is a bug, and the read-only env cannot install.

## 10. Environment probe — `scatac_env_probe`

**Call.** `scatac_env_probe(module=<1-6>, write_panel=<project>/panel.json)`.
Call it before promising motif enrichment, gene activity, integration or any
optional package. `write_panel` writes a starter marker-panel JSON for §6.

**Read.** `environment`; `packages_total`/`packages_installed`; `packages[]`
(`{name, version, installed}`); `modules[]` (`{module, missing_required,
missing_optional, core_ready}`); `core_ready_modules`; `blocked_modules`;
`panel_written`; `notes`.

**Decide.**
- A `core_ready` module can run here; a `blocked` module names the missing
  package and the project-environment next step. The core chain (snapatac2,
  scanpy, anndata, polars, pyarrow, macs3, leidenalg, igraph, umap-learn,
  scikit-learn, scipy, numpy, pandas, h5py, plotly, matplotlib) is shipped.
- chromVAR / ArchR / Signac (R), episcanpy (numpy 2 incompatible), scvi-tools /
  magic-impute (pip-only) are gaps: report them as facts with an `env_request`
  project-environment next step, never an in-place install.

**Red lines.** The plugin environment is read-only and prebuilt — never install
into it. A gap is a fact plus a next step, not a silent skip and not an
impersonated capability.
*Why:* managed environments are reproducible and immutable; optional heavy tools
belong in a user project environment.

## 11. Red-line checklist

Before delivering, walk the table in
[`red-line-cheatsheet.md`](red-line-cheatsheet.md) (12 rows). The short version:

1. Raw accessibility first — peaks/diff/scrublet on fragments/counts, never
   corrected/scaled/integrated values.
2. Every backed AnnData is closed before exit (script saves+closes; notebook
   ends with a close cell).
3. The analysis unit is the cell (barcode), not a spot/bin; tiles/peaks are
   genome features, not units.
4. Doublets inflate accessibility — scrublet matters; report the doublet rate.
5. TSSe is the ATAC signal-to-noise metric; thresholds come from the
   distribution; `kept_fraction` < 0.5 → stop.
6. Tile bin / peak width / per-cluster calling are biological decisions — state
   them.
7. Gene activity is a proxy, not RNA — annotation/marker viz only, no DE claims.
8. Label transfer needs a scRNA reference (from the single-cell plugin); without
   one, only gene-activity inference, labelled as such.
9. Motif enrichment is sequence prediction, needs FASTA + motif DB; chromVAR is
   an R gap, not shipped — do not impersonate.
10. The sample is the experimental unit; ≥2 biological replicates per condition
    for a condition-level claim.
11. Integration is for clustering/viz only; do not over-correct; corrected values
    never enter diff.
12. No alignment/Cell Ranger in the plugin; missing genome resource/package →
    gap report + project-env next step, never an in-place install.

A step carrying a "skipped" note did not run — never read its absence as a clean
result. Quote the tool's JSON numbers; do not re-derive them.

## Sources and licenses

Knowledge is an original distillation of the following open sources; no upstream
text, figures or data are copied.

| Source | URL | License |
| --- | --- | --- |
| SnapATAC2 documentation & tutorials (Zhang et al.) | https://github.com/kaizhang/SnapATAC2 | MIT |
| Single-cell best practices (chromatin accessibility chapter) | https://www.sc-best-practices.org/ | Apache 2.0; Heumos et al., Nat Rev Genet 2023 |
| ENCODE Blacklist (Amemiya, Kundaje, Boyle) | doi:10.1038/s41598-019-45839-z | Sci Rep 9, 9354 (2019); Boyle-Lab/Blacklist ENCODE v2 lists |

The ENCODE hg38/mm10 blacklist regions are inlined (base64-gzip) in
`assets/fragments/blacklist.py` so both the script tools and the rendered
notebook filter artifact regions offline. snapatac2 bundles chrom_sizes for
hg38/mm10/hg19/mm39 offline; gene annotation (GFF3) and genome FASTA are lazy
downloads provided by the user via `input-path` or downloaded once into a project
directory, never into the read-only plugin environment. The code is original to
this plugin and runs in the private `plugin:scatac` environment.