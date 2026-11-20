# Single-cell analysis playbook (scanpy route)

The agent's operating manual. Each stage below is a decision procedure:
**Call** the tool, **Read** the named JSON fields and figures, **Decide** the
parameters from what you read, and respect the **Red lines**. The one-line
*why* gives the principle; the hard rules (route choice, renderer/agent split,
one reading for both paths, lognorm snapshot) live in SKILL.md `Rules`; the
learner-facing short version lives in `references/step-notes/` (one note per
step). Sources
and licenses are at the end.

Every number here is a **starting point for 10x-scale data, not a default to
accept**. Read the distributions first, adapt, and say in your answer which
value you chose and why. Quote the tool's JSON numbers instead of re-deriving
them.

## 0. The decision loop

Chain: `sc_py_inspect` → `sc_py_qc` → `sc_py_preprocess` →
(`sc_py_integrate` when batches separate) → `sc_py_markers` →
`sc_py_annotate` → `sc_py_pseudobulk`; `sc_py_trajectory` when the question
is a process, not groups. `sc_py_notebook` is the parallel teaching path;
feed it the *same* threshold readings (flags or `config` JSON). Skip a step
only when `sc_py_inspect`'s `candidate_steps` says its result already exists.

Loop at every stage: call with loose/starting values → read the JSON →
decide → re-call with the decided values → verify `after`, `warnings`,
`notes`. A default that reached the output untouched is a decision you did
not make.

Call `sc_py_env_probe` before promising scVI, CellTypist, enrichment or
spatial work: the plugin environment is read-only; a gap is reported as a
fact with a project-environment next step, never an in-place install.

## 1. Inspect the input — `sc_py_inspect`

**Call.** `sc_py_inspect(input, output=<project>/loaded.h5ad)`. Keep
`min_features=0` here: gating belongs to `sc_py_qc`, which reports per-gate
removals; filtering at load time deletes cells silently. For 10x input choose
`var_names` (`gene_symbols` default; `gene_ids` when the downstream panel is
Ensembl-based) — the identifier type decides every later regex and panel.

**Read.** `shape`; `counts_present` (false → re-import the raw matrix before
QC, seurat_v3 HVG or pseudobulk); `value_class`; `metadata_columns` (choose
`sample_column`/`batch_key`: the donor/library column, not a cluster column);
`candidate_steps`; `notes`.

**Decide.** `sample_column` and `batch_key` for everything downstream; which
steps to skip.

**Red lines.** No QC, counts-flavor HVG or pseudobulk without a `counts`
layer.
*Why:* raw counts are the only comparable evidence; normalized values cannot
be re-normalized or summed into a counts model.

## 2. Quality control — `sc_py_qc`

**Call 1 (reading).** Leave `feature_low/high` and `count_low/high` unset
(a None bound skips that gate) and pass a loose `mt_max` (e.g. 100) — the
adapter defaults `mt_max` to 10, which would already gate the distribution
you are trying to read. Keep `detect_doublets` off for a pure reading call:
the filter step always runs, and `true` removes Scrublet-flagged cells even
with every numeric bound unset (None bounds skip only their own gates).
Decide doublet removal for the gating run. Set `figures_dir` and look at the
before-violins **per sample** (`before.per_sample`), not one pooled plot.
`detect_doublets=true` runs Scrublet (shipped in this environment).

**Read.** `distribution[]` for `n_genes_by_counts`, `total_counts`,
`pct_counts_mt`: `min/q01/median/q99/max`. `doublets` (`method=scrublet`,
`expected_rate`, `singlet/doublet`, or `error` — Scrublet can fail on tiny or
very sparse inputs; the report says so instead of pretending zero doublets).

**Decide.**
- `feature_low`: the valley between the empty-droplet shoulder and the cell
  peak; start ~200 genes for 10x, lower for shallow libraries. Cross-check
  with an adaptive rule (OSCA `isOutlier`: 3 MADs below the median on
  log-scale) and confirm the gate does not cut the cell peak.
- `count_low`: same logic on depth; start a few hundred UMIs, adapt to the
  `total_counts` distribution.
- `feature_high` / `count_high`: the upper (doublet) shoulder — set from
  `q99`/`max` and the violin, not by default.
- `mt_max`: tissue-dependent starting points — cell lines/PBMC 5–10, solid
  tumour 10–20, heart/kidney/muscle 20–50. Read the `pct_counts_mt`
  distribution (damage is a shoulder, not a round number) before fixing the
  cap.
- `mt_pattern`: `^MT-` human, `^mt-` mouse. If `pct_counts_mt` median is 0,
  the regex or identifier type is wrong — fix that first; a gate that never
  fires is not QC.
- `doublet_expected_rate`: a loading estimate, ~0.8% per 1000 cells loaded
  (the adapter's default), not a free knob. Run per sample when libraries
  were loaded separately. Heterotypic doublets invent "hybrid" clusters.

**Call 2 (gating).** Pass the decided values. **Read** `filters[]`
(`filter`, `removed` — trace a surprising total to one gate),
`after.kept_fraction`, `warnings`, and the after-violins. Scoring alone
removes nothing; the filter step applies the doublet call.

**Red lines.** `kept_fraction` < 0.5 → stop and re-read the distributions;
zero cells kept → relax. Libraries that differ in depth/stress get per-sample
gates. QC runs on raw counts only.
*Why:* gene count, depth and mitochondrial fraction carry distinct signals
(empty droplet / debris / doublet); adaptive gates beat round numbers.

## 3. Normalize, reduce, cluster — `sc_py_preprocess`

**Read.** `top_variable_features` (ribosomal/mitochondrial genes or one
population monopolizing the list = warning); `variance`
(`pcs`, `variance_explained[]`, `first_ten_percent`) plus the `pca-elbow`
figure; `clusters[]` per resolution column (`clusters`, `sizes[]`);
`dimplot-clusters`; `notes`. The step snapshots pre-scale log-normalized
expression to `layers["lognorm"]` — the expression of record for markers and
panel scoring.

**Decide.**
- `hvg_flavor`: `seurat_v3` default, on the **counts** layer (needs
  scikit-misc, shipped). A counts flavor on log values silently mis-ranks.
- `n_variable_features`: start 2000; raise toward 3000 when rare populations
  matter. With several samples set `batch_key` so HVGs are ranked per batch —
  otherwise one sample chooses every gene.
- `n_pcs` 30 to compute; `dims_use` from the elbow **plus a few PCs beyond the
  bend** so rare populations survive (start 10–30, read `variance_explained`).
- `resolution`: sweep `[0.4, 0.8, 1.2]`, compare `sizes` and stability across
  the sweep; pick the column whose clusters are stable and sensibly sized,
  and use that column as `group_by` downstream.
- `regress_vars`: off by default; regressing covariates can remove real
  biology with them. `normalization` here is lognormalize only — SCT lives on
  the Seurat route.

**Red lines.** A cluster that is one donor is a batch artifact until proven
otherwise (§4). Tiny clusters (tens of cells) → suspect doublets or
over-splitting before biology. Never cluster on the 2D embedding; subcluster
with fresh HVG/PCA on the subset.
*Why:* HVGs rank variance beyond the mean–variance trend; PCA orders signal
before noise; resolution is a granularity dial, clusters are hypotheses, not
cell types.

## 4. Batch integration — `sc_py_integrate`

**When.** More than one sample and the batch-colored UMAP from
`sc_py_preprocess` shows donors separating. One sample → skip.

**Call.** `method=harmony` default (corrects the PCA embedding, stored as
`obsm["X_pca_integrated"]`); `combat` when you want expression-level
correction (runs on the pre-scale lognorm snapshot, then recomputes PCA);
`bbknn` is a documented gap (no python 3.13 build) — the tool reports it, use
a project environment instead. `batch_key` = the donor/library column.

**Read.** `figures` `batch-before`/`batch-after`; `clusters[]` (re-clustered
on the corrected embedding); `notes` (which space was corrected); `method`,
`batches`.

**Decide.** Judge on the before/after pair: samples must mix **within each
cell type** while markers and known populations survive (scIB weights batch
removal ~40% / biology conservation ~60% — mixing alone is not success). A
whole biological group vanishing after correction = over-correction; try the
other method or weaken the correction. Downstream markers/annotation read
the corrected clusters.

**Red lines.** Corrected values NEVER feed differential expression — test on
uncorrected values with sample as the unit (§6). Without a batch column the
chain continues uncorrected and says so; do not invent a batch from cell-type
labels.
*Why:* methods differ in what they correct (embedding vs expression vs
graph); the verdict is visual plus biology conservation, not a score.

## 5. Marker genes — `sc_py_markers`

**Call.** `group_by` = the resolution column chosen in §3/§4; `test=wilcoxon`,
`logfc_threshold=0.25`, `min_pct=0.1`, `only_positive=true` as starting
points; `tables_dir` + `figures_dir` set.

**Read.** `groups[]` (top genes with log2 fold changes), `total_markers`, the
per-group CSVs in `tables` (they carry the detection fractions `pts` and
adjusted p-values — the report's top-10 is a summary, the CSV is the
evidence), the heatmap, `notes` — including whether the test ran on the
`lognorm` snapshot ("no lognorm snapshot; markers tested on the current X"
is a warning to fix, not to ignore). `total_markers=0` → lower the thresholds
or the clusters do not differ.

**Decide.** A marker is a gene with high detection in the cluster **and** low
detection in the rest (specific), not merely a large fold change. Read several
candidates per cluster and check the heatmap for cross-cluster exclusivity
before naming anything. Unequal cluster sizes distort rank tests — downsample
the large cluster when comparing.

**Red lines.** Adjusted p-values are descriptive only: clusters were defined
from the same data (double dipping), never report them as inferential
statistics. Tests never run on scaled or batch-corrected values.
`avg_log2FC`/logfoldchanges are this comparison's fold change, not a portable
effect size.
*Why:* Wilcoxon is the best-calibrated rank test here; detection fractions
say how specific a gene is, significance alone does not.

## 6. Cell-type annotation — `sc_py_annotate`

**Call.** Build a panel JSON (cell type → gene symbols; starter via
`sc_py_env_probe write_panel`), coarse lineages first, subtypes second.
`method=panel`; `mapping` JSON only when the user has ground truth.

**Read.** `scores[]` per cluster (per-panel scores, `genes_missing`;
`best=null` means below `min_score`), `assignments`, `notes` (unassigned
cluster list; absent-panel-gene warning), the dot plot.

**Decide.** `min_score` starts 0; look at the winner–runner-up margin in
`scores` — a narrow win between two similar panels is not an identification,
raise `min_score` or leave the cluster unassigned. Panel symbols must match
the object's identifier type and species case.

**Red lines.** Zero usable panel genes = identifier mismatch, not a biology
result. Report how many clusters got each label and how many stayed
unassigned — unassigned is information, not a failure to hide. `celltypist`
needs the package and a model, absent from the plugin environment: report the
gap (and it expects log1p 1e4 input). Never propagate a coarse label into
subtype claims.
*Why:* panel-mean-minus-rest keeps panel sizes comparable; a name must rest
on stated evidence.

## 7. Pseudobulk handoff — `sc_py_pseudobulk`

**Call.** `sample_column` = the biological replicate (donor), `group_column`
= the condition; `min_cells=10`, `min_cells_expressing=3` as starting points.

**Read.** `design` (one row per sample×group unit with its cell count),
`units`, `dropped_genes`, `notes` (thin groups, skipped units), and
`wrapper_hint` (`counts`, `samplesheet`, `sample_id_col`, `gene_id_col`,
`contrast_variable`).

**Decide.** Pass `wrapper_hint` verbatim to the repository's DESeq2/edgeR/
limma wrapper — the layout already matches, and it is the only sanctioned
cross-route handoff besides `references/cross-route-handoff.md`. Power comes from more
donors, not more cells per donor.

**Red lines.** Raw counts only (summing normalized/scaled values invalidates
the counts model). Fewer than `min_samples_per_group` (default 2) replicates
per group → no condition-level claim; say so. Cell-level tests on one donor's
cells are pseudoreplication.
*Why:* the donor, not the cell, is the experimental unit.

## 8. Trajectory — `sc_py_trajectory`

**When.** The question is a process (differentiation, activation), not groups.
Needs a neighbor graph: run after `sc_py_preprocess` or `sc_py_integrate`.

**Call.** `root_cluster` = the cluster that **starts** the process, chosen
from prior biology (stem/progenitor markers) — the default (first level) is a
decision you did not make; name your root and why. `n_dcs=10` starting point.

**Read.** `root_cluster`, `pseudotime` (`min`, `max`, `unreachable_cells`),
the PAGA graph and pseudotime UMAP figures, `notes`.

**Decide.** Read PAGA edges as hypotheses about transitions and validate
against known marker gradients along pseudotime. `unreachable_cells` (inf
pseudotime) are disconnected components, not "late" cells — report them.
Pseudotime is a relative rank within one path: not a clock, not comparable
across paths or datasets.

**Red lines.** Parameters can manufacture the trajectory you expected;
validation against independent biology is mandatory. A trajectory through
over-corrected integrated data can be an artifact — check §4 first.
*Why:* PAGA abstracts the neighbour graph to cluster connectivity; DPT ranks
cells by diffusion distance from the root.

## 9. The notebook path — `sc_py_notebook`

**Decide.** Preset (`module` 1 QC, 2 +clustering, 3 +integration, 4 +markers,
5 +pseudobulk, 6 +annotation, 7 +trajectory) or `steps` with automatic
dependency closure. Write the §2–§4 readings into `config` JSON
(`min_genes`, `max_counts`, `mt_max`, `mt_pattern`, `n_hvg`, `n_pcs`, `dims`,
`resolution`, `seed`, `batch_key`, `batch_method`, `hvg_flavor`) so the
notebook and the script chain agree by construction. `self_check=true` before
delivery. `inject`/`merge` sources must match the notebook language.

**Read.** The report's `steps`, `code_cells`, `self_check`; the notebook's
`phi_render` metadata (parameters with sources, fragment sha256, data
fingerprint) is the audit record.

**Red lines.** The renderer fills; you decide. Every config value carries a
source tag — a `[default]` tag on a threshold means you skipped a decision.

## 10. Red-line checklist

- No counts layer → re-import raw before QC/seurat_v3/pseudobulk.
- `kept_fraction` < 0.5 or zero cells → stop, re-read distributions.
- `pct_counts_mt` median 0 → wrong regex/identifier type, fix before gating.
- Markers and panel scoring read `layers["lognorm"]`; never scaled or
  corrected values; "no lognorm snapshot" note → fix the chain.
- Corrected embeddings cluster; corrected values never test.
- Marker p-values descriptive; fold changes not portable.
- Unassigned clusters and absent panel genes are reported, not hidden.
- <2 replicates per group → no condition claim.
- Trajectory root named from biology; inf pseudotime reported.
- Gaps (bbknn, celltypist, scvi-tools, gseapy, squidpy…) → project
  environment, never an in-place install.

## Sources and licenses

| Source | URL | License |
| --- | --- | --- |
| Orchestrating Single-Cell Analysis (OSCA) | https://bioconductor.org/books/release/OSCA/ | CC BY 4.0 |
| Single-cell best practices | https://www.sc-best-practices.org/ | Apache 2.0; cite Heumos, Schaar, Lance et al., Nat Rev Genet 2023, DOI 10.1038/s41576-023-00586-w |
| Cambridge single-cell course (Hemberg lab) | https://www.singlecellcourse.org/ | GPL-3 (course material; not OSCA) |
| NBIS workshop-scRNAseq | https://nbisweden.github.io/workshop-scRNAseq/ | CC-BY-4.0 per repository LICENSE (site footer says GPL-3; cite the repo) |
| scNotebooks curriculum | https://github.com/integrativebioinformatics/scNotebooks | CC BY 4.0; Rojas-Hidalgo et al., Nat Genet 2026 |

All prose is original, distilled from these open textbooks; no upstream text
is copied. Topic chapters: best-practices preprocessing (qc, doublet_detection,
normalization, hvg, dimensionality_reduction, clustering, integration),
annotation (marker_identification, cell_type_annotation), downstream
(differential_expression, trajectory_inference); OSCA quality-control,
normalization, clustering, marker-detection, cell-type-annotation,
multi-sample-comparisons.