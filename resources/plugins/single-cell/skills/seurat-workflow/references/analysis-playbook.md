# Single-cell analysis playbook (Seurat route)

The agent's operating manual. Each stage below is a decision procedure:
**Call** the tool, **Read** the named JSON fields and figures, **Decide** the
parameters from what you read, and respect the **Red lines**. The one-line
*why* gives the principle; the hard rules (route choice, renderer/agent split,
one reading for both paths) live in SKILL.md `Rules`; the learner-facing short
version lives in `references/step-notes/` (one note per step). Sources and
licenses are at the
end.

Every number here is a **starting point for 10x-scale data, not a default to
accept**. Read the distributions first, adapt, and say in your answer which
value you chose and why. Quote the tool's JSON numbers instead of re-deriving
them.

## 0. The decision loop

Chain: `sc_inspect` → `sc_qc` → `sc_preprocess` → `sc_markers` →
`sc_annotate` → `sc_pseudobulk`. `sc_notebook` is the parallel teaching path;
feed it the *same* threshold readings (flags or `config` JSON). Skip a step
only when `sc_inspect`'s `candidate_steps` says its result already exists.

Loop at every stage: call with loose/starting values → read the JSON →
decide → re-call with the decided values → verify `after`, `warnings`,
`notes`. A default that reached the output untouched is a decision you did
not make.

Call `sc_env_probe` before promising trajectory, cell–cell communication,
spatial or scATAC work: the plugin environment is read-only; a gap is reported
as a fact with a project-environment next step, never an in-place install.

## 1. Inspect the input — `sc_inspect`

**Call.** `sc_inspect(input, output=<project>/loaded.rds)`. Keep
`min_features=0` here: gating belongs to `sc_qc`, which reports per-gate
removals; `min.features` at load time deletes cells silently.

**Read.** `kind`, `shape`; `counts_present` (false → re-import the raw matrix
before QC or pseudobulk); `value_class`; `metadata_columns` (choose the
`sample_column`: the donor/library column, not a cluster column);
`candidate_steps` (which of qc/markers/annotate/pseudobulk the object is ready
for); `notes`.

**Decide.** `sample_column` for everything downstream; which steps to skip.

**Red lines.** No QC or pseudobulk on an object without a `counts` layer.
*Why:* raw counts are the only comparable evidence; normalized values cannot
be re-normalized or summed into a counts model.

## 2. Quality control — `sc_qc`

**Call 1 (reading).** Leave `feature_low/high` and `count_low/high` unset
(a NULL bound skips that gate) and pass a loose `mt_max` (e.g. 100) — the
adapter defaults `mt_max` to 10, which would already gate the distribution
you are trying to read. Keep `detect_doublets` off for a pure reading call:
the filter step always runs, and `true` removes flagged doublets even with
every numeric bound unset (NULL bounds skip only their own gates). The
scDblFinder gap comes from `sc_env_probe`; `doublets.error` appears once
detection is attempted in the gating run. Set `figures_dir` and look at the
before-violins **per sample** (`before.per_sample`), not one pooled plot.

**Read.** `distribution[]` for `nFeature_RNA`, `nCount_RNA`, `percent.mt`:
`min/q01/median/q99/max`. `doublets` (`method`, `expected_rate`,
`singlet/doublet`, or `error`). The plugin environment lacks scDblFinder:
`doublets.error` says so — report the gap, do not install.

**Decide.**
- `feature_low`: the valley between the empty-droplet shoulder and the cell
  peak; start ~200 genes for 10x, lower for shallow libraries. Cross-check
  with an adaptive rule (OSCA `isOutlier`: 3 MADs below the median on
  log-scale) and confirm the gate does not cut the cell peak.
- `count_low`: same logic on depth; start a few hundred UMIs, adapt to the
  `nCount_RNA` distribution.
- `feature_high` / `count_high`: the upper (doublet) shoulder — set from
  `q99`/`max` and the violin, not by default.
- `mt_max`: tissue-dependent starting points — cell lines/PBMC 5–10, solid
  tumour 10–20, heart/kidney/muscle 20–50. Read the `percent.mt` distribution
  (damage is a shoulder, not a round number) before fixing the cap.
- `mt_pattern`: `^MT-` human, `^mt-` mouse. If `percent.mt` median is 0, the
  regex or identifier type is wrong (Ensembl IDs need a different pattern) —
  fix that first; a gate that never fires is not QC.
- `doublet_expected_rate`: a loading estimate, ~0.8% per 1000 cells loaded
  (the adapter's default), not a free knob. Heterotypic doublets invent
  "hybrid" clusters; homotypic ones are mostly harmless.

**Call 2 (gating).** Pass the decided values. **Read** `filters[]`
(`filter`, `removed` — trace a surprising total to one gate),
`after.kept_fraction`, `warnings`, and the after-violins.

**Red lines.** `kept_fraction` < 0.5 → stop and re-read the distributions;
zero cells kept → relax. Libraries that differ in depth/stress get per-sample
gates. QC runs on raw counts only.
*Why:* gene count, depth and mitochondrial fraction carry distinct signals
(empty droplet / debris / doublet); adaptive gates beat round numbers.

## 3. Normalize, reduce, cluster — `sc_preprocess`

**Read.** `top_variable_features` (ribosomal/mitochondrial genes or one
population monopolizing the list = warning); `variance`
(`pcs`, `variance_explained[]`, `first_ten_percent`) plus the `pca-elbow`
figure; `clusters[]` per resolution column (`clusters`, `sizes[]`);
`dimplot-clusters`; `notes`.

**Decide.**
- `normalization`: `lognormalize` default; `sct` when depth varies wildly or
  the object feeds integration (SCT lives on this route).
- `n_variable_features`: start 2000; raise toward 3000 when rare populations
  matter. Too few merges rare types; too many re-admits noise.
- `n_pcs` 30 to compute; `dims_use` from the elbow **plus a few PCs beyond the
  bend** so rare populations survive (start 10–30, read `variance_explained`).
- `resolution`: sweep `[0.4, 0.8, 1.2]`, compare `sizes` and stability across
  the sweep (clustree/ARI if you render them); pick the column whose clusters
  are stable and sensibly sized, and use that column as `group_by` downstream.
- `regress_vars`: off by default; regressing covariates can remove real
  biology with them.

**Red lines.** A cluster that is one donor is a batch artifact until proven
otherwise (§7). Tiny clusters (tens of cells) → suspect doublets or
over-splitting before biology. Never cluster on the 2D embedding; subcluster
with fresh HVG/PCA on the subset.
*Why:* HVGs rank variance beyond the mean–variance trend; PCA orders signal
before noise; resolution is a granularity dial, clusters are hypotheses, not
cell types.

## 4. Marker genes — `sc_markers`

**Call.** `group_by` = the resolution column chosen in §3; `test=wilcox`,
`logfc_threshold=0.25`, `min_pct=0.1`, `only_positive=true` as starting
points; `tables_dir` + `figures_dir` set.

**Read.** `groups[]` (`top_genes`, `top_log2fc`), `total_markers`, the
per-group CSVs in `tables` (they carry `pct.1`, `pct.2`, `p_val_adj` — the
report's top-10 is a summary, the CSV is the evidence), the heatmap, `notes`
(`total_markers=0` → lower the thresholds or the clusters do not differ).

**Decide.** A marker is a gene with high `pct.1` **and** low `pct.2`
(specific), not merely a large fold change. Read several candidates per
cluster and check the heatmap for cross-cluster exclusivity before naming
anything. Unequal cluster sizes distort rank tests — downsample the large
cluster when comparing.

**Red lines.** `p_val_adj` is descriptive only: clusters were defined from the
same data (double dipping), never report it as inferential statistics. Tests
never run on scaled or batch-corrected values. `avg_log2FC` is this
comparison's fold change, not a portable effect size.
*Why:* Wilcoxon is the best-calibrated rank test here; detection fractions
say how specific a gene is, significance alone does not.

## 5. Cell-type annotation — `sc_annotate`

**Call.** Build a panel JSON (cell type → gene symbols; starter via
`sc_env_probe write_panel`), coarse lineages first, subtypes second.
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
unassigned — unassigned is information, not a failure to hide. `singler`
needs SingleR/celldex, absent from the plugin environment: report the gap.
Never propagate a coarse label into subtype claims.
*Why:* panel-mean-minus-rest keeps panel sizes comparable; a name must rest
on stated evidence.

## 6. Pseudobulk handoff — `sc_pseudobulk`

**Call.** `sample_column` = the biological replicate (donor), `group_column`
= the condition; `min_cells=10`, `min_cells_expressing=3` as starting points.

**Read.** `design` (one row per sample×group unit with its cell count),
`units`, `dropped_genes`, `notes` (thin groups, skipped units), and
`wrapper_hint` (`counts`, `samplesheet`, `sample_id_col`, `gene_id_col`,
`contrast_variable`).

**Decide.** Pass `wrapper_hint` verbatim to the repository's DESeq2/edgeR/
limma wrapper — the layout already matches. Power comes from more donors,
not more cells per donor.

**Red lines.** Raw counts only (summing normalized/scaled values invalidates
the counts model). Fewer than `min_samples_per_group` (default 2) replicates
per group → no condition-level claim; say so. Cell-level tests on one donor's
cells are pseudoreplication.
*Why:* the donor, not the cell, is the experimental unit.

## 7. Batch effects — capability boundary of this route

Diagnose with the notebook's `batch-compare` step (UMAP colored by
`sample_column`): one donor sitting in one corner is the uncorrected-batch
signature. **The R script chain has no integration tool**, and the plugin
environment lacks `harmony`. Options, in order: hand the object to the Python
route (`sc_py_integrate`, see scanpy-workflow `references/cross-route-handoff.md`), or
use a project environment with harmony/Seurat integration. Never install
inside the plugin environment.

**Red lines.** Corrected values never feed marker tests or differential
expression — test on uncorrected values with sample as the unit. Judge
correction by before/after embeddings: mixing **and** biology conservation;
a vanished biological group means over-correction.

## 8. The notebook path — `sc_notebook`

**Decide.** Preset (`module` 4 QC→pseudobulk, 5 +batch compare, 8 +ADT,
9 +clonotypes) or `steps` with automatic dependency closure. Write the §2–§3
readings into `config` JSON (`min_genes`, `max_counts`, `mt_max`,
`mt_pattern`, `n_hvg`, `n_pcs`, `dims`, `resolution`, `seed`) so the notebook
and the script chain agree by construction. `self_check=true` before
delivery. `inject`/`merge` sources must match the notebook language.

**Read.** The report's `steps`, `code_cells`, `self_check`; the notebook's
`phi_render` metadata (parameters with sources, fragment sha256, data
fingerprint) is the audit record.

**Red lines.** The renderer fills; you decide. Every config value carries a
source tag — a `[default]` tag on a threshold means you skipped a decision.

## 9. Red-line checklist

- No counts layer → re-import raw before QC/pseudobulk/HVG-on-counts.
- `kept_fraction` < 0.5 or zero cells → stop, re-read distributions.
- `percent.mt` median 0 → wrong regex/identifier type, fix before gating.
- Never test or sum scaled/corrected values; markers and DE read pre-scale
  log-normalized expression (this route: the `RNA` assay data slot).
- Marker p-values descriptive; `avg_log2FC` not portable.
- Unassigned clusters and absent panel genes are reported, not hidden.
- <2 replicates per group → no condition claim.
- Gaps (scDblFinder, harmony, SingleR, monocle3, liana…) → project
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