# Spatial transcriptomics analysis playbook (Python route)

The agent's operating manual. Each stage below is a decision procedure:
**Call** the tool, **Read** the named JSON fields and figures, **Decide** the
parameters from what you read, and respect the **Red lines**. The one-line
*why* gives the principle; the hard rules (platform/bin first, renderer/agent
split, one reading for both paths, counts/lognorm layers) live in SKILL.md
`Rules`; the learner-facing short version lives in `references/step-notes/`
(one note per step). Sources and licenses are at the end.

Every number here is a **starting point for Visium/Stereo-seq data, not a
default to accept**. Read the distributions first, adapt, and say in your answer
which value you chose and why. Quote the tool's JSON numbers instead of
re-deriving them.

## 0. The decision loop

Chain: `st_py_inspect` → `st_py_qc` → `st_py_preprocess` → then branch on the
question:

- **What are the domains?** → `st_py_markers` → `st_py_annotate`.
- **Which genes carry spatial structure / which domains are neighbours?** →
  `st_py_spatial` (SVG + neighbourhood enrichment/co-occurrence).
- **What cell types make up each spot/bin?** → `st_py_deconvolve` (needs a
  single-cell reference).
- **Which domains might signal to each other?** → `st_py_communication`.
- **Does spatial organisation differ between conditions?** → `st_py_spatial`
  with `group_column` set (cross-sample differential).

`st_py_notebook` is the parallel teaching path; feed it the *same* threshold
readings (flags or `config` JSON). Skip a step only when `st_py_inspect`'s
`candidate_steps` says its result already exists.

Loop at every stage: call with loose/starting values → read the JSON → decide →
re-call with the decided values → verify `after`, `warnings`, `notes`. A default
that reached the output untouched is a decision you did not make.

Call `st_py_env_probe` before promising cell2location, tangram, stereopy, liana
or any model-based tool: the plugin environment is read-only; a gap is reported
as a fact with a project-environment next step, never an in-place install.

## 1. Inspect the input — `st_py_inspect`

**Call.** `st_py_inspect(input, output=<project>/loaded.h5ad)`. Keep
`min_genes=0`/`min_units=0` here: gating belongs to `st_py_qc`, which reports
per-gate removals; filtering at load time deletes units silently. For a count
table choose `var_names` (`gene_symbols` default; `gene_ids` when the downstream
panel is Ensembl-based) — the identifier type decides every later regex and
panel. For Stereo-seq pass `bin_size` or `cellbin` deliberately (see Decide).

**Read.** `platform` (visium / visium-hd / stereoseq, and `bin_auto`/`bin_size`/
`cellbin` provenance); `shape` (`units`, `genes`, `layers`, `obsm`, `obsp`,
`snapshot`); `counts_present` (false → re-import the raw matrix before QC,
seurat_v3 HVG or deconvolution); `coordinates_present` (false → only non-spatial
analysis is possible); `value_class`; `metadata`/`metadata_columns` (choose
`sample_column`/`group_column`: the donor/condition column, not a cluster
column); `candidate_steps`; `notes`.

**Decide.**
- **Platform and analysis unit first.** Visium → spot (several cells per spot).
  Visium HD → square bin (8/16µm) or seg (cells); inspect says which. Stereo-seq
  → square bin or cellbin. State the unit you analyse.
- **Stereo-seq bin size is a biological decision.** bin20≈10µm, bin50≈25µm;
  larger bins for tissue domains, cellbin or small bins for cell-level questions.
  A `bin_auto` note means the reader picked one — override it explicitly and say
  why. cellbin segmentation is done upstream by SAW; the plugin only reads it.
- `sample_column` and `group_column` for everything downstream; which steps to
  skip.

**Red lines.** No QC, counts-flavor HVG or deconvolution without a `counts`
layer. No SVG/neighbourhood/communication/spatial maps without
`obsm["spatial"]` coordinates and a spatial adjacency graph.
*Why:* raw counts are the only comparable evidence; spatial statistics are
defined on tissue coordinates, not on expression space.

## 2. Quality control — `st_py_qc`

**Call 1 (reading).** Leave `min_genes`/`min_counts`/`max_counts` unset (a None
bound skips that gate) and pass a loose `mt_max` (e.g. 100) so you read the
distribution before gating it. Set `figures_dir` and look at the before-violins
and the **spatial** QC maps (`qc-spatial-before`): low-quality units often sit
at tissue edges or in blank regions, which a pooled violin hides.

**Read.** `distribution[]` for `n_genes_by_counts`, `total_counts`,
`pct_counts_mt`: `min/q01/median/q99/max`. `platform`/`bin`. `before.per_sample`
when a sample column exists.

**Decide.**
- `min_genes`: spatial units are shallower than dissociated cells — start ~100
  for Visium spots, lower for small Stereo-seq bins; read the valley between the
  empty-unit shoulder and the tissue peak. Cross-check with an adaptive rule
  (3 MADs below the median log-depth) and confirm the gate does not cut the peak.
- `min_counts`/`max_counts`: from the `total_counts` distribution; the upper
  bound catches multi-cell aggregates only when the unit is meant to be smaller.
- `mt_max`: tissue-dependent — spatial tissue tolerates more than dissociated
  cells; start 20 and read the `pct_counts_mt` shoulder before fixing the cap.
- `mt_pattern`: `^MT-` human, `^mt-` mouse. If `pct_counts_mt` is absent or its
  median is 0, the regex or identifier type is wrong — fix that first; a gate
  that never fires is not QC (the report warns when no gene matched).

**Call 2 (gating).** Pass the decided values. **Read** `filters[]` (`filter`,
`removed` — trace a surprising total to one gate), `after.kept_fraction`,
`warnings`, and the after-violins/spatial maps.

**Red lines.** `kept_fraction` < 0.5 → stop and re-read the distributions (the
bin may be too small); zero units kept → relax. Samples that differ in
depth/stress get per-sample gates. QC runs on raw counts only.
*Why:* gene count, depth and mitochondrial fraction carry distinct signals;
spatial QC adds tissue coverage — a unit can pass numerically yet lie off-tissue.

## 3. Normalize, reduce, cluster, build the spatial graph — `st_py_preprocess`

**Read.** `top_variable_features` (ribosomal/mitochondrial genes or one
population monopolizing the list = warning); `variance` (`pcs`,
`variance_explained[]`, `first_ten_percent`) plus the `pca-elbow` figure;
`clusters[]` per resolution column (`clusters`, `sizes[]`); `spatial_graph`
(`method`, `edges`, `mean_degree`); `dimplot-clusters` and `spatial-clusters`;
`notes`. The step snapshots pre-scale log-normalized expression to
`layers["lognorm"]` — the expression of record for markers, panel scoring, SVG
and communication.

**Decide.**
- `hvg_flavor`: `seurat_v3` default, on the **counts** layer (needs scikit-misc,
  shipped). A counts flavor on log values silently mis-ranks.
- `n_hvg`: start 2000; raise toward 3000 when rare domains matter.
- `n_pcs` 30 to compute; `dims` from the elbow **plus a few PCs beyond the bend**
  so rare domains survive (read `variance_explained`).
- `resolution`: sweep `[0.4, 0.8, 1.2]` (repeat the flag), compare `sizes` and
  stability; pick the column whose clusters are stable and sensibly sized, and
  use that column as `group_by` downstream.
- `neighbor_method` + `n_neighs`/`radius`: the spatial adjacency geometry.
  Visium hex lattice → `grid` (n_neighs=6) or `delaunay`; Stereo-seq bins or
  irregular layouts → `knn` or `radius`. Read `spatial_graph.mean_degree`: too
  sparse misses neighbours, too dense bridges gaps. State the method — it changes
  every downstream spatial conclusion.

**Red lines.** A cluster that is one sample is a batch artifact until proven
otherwise. Tiny clusters (tens of units) → suspect over-splitting or empty
regions before biology. Never cluster on the 2D embedding. The spatial graph is
built on coordinates, NOT on the expression PCA.
*Why:* HVGs rank variance beyond the mean–variance trend; PCA orders signal
before noise; resolution is a granularity dial; the adjacency graph is the
substrate of all spatial statistics.

## 4. Marker genes — `st_py_markers`

**Call.** `group_by` = the resolution column chosen in §3; `test_use=wilcoxon`,
`logfc_threshold=0.25`, `min_pct=0.1`, `only_positive=true` as starting points;
`tables_dir` + `figures_dir` set.

**Read.** `groups[]` (top genes with log2 fold changes), `total_markers`, the
per-group CSVs in `tables` (they carry detection fractions and adjusted p-values
— the report's top-10 is a summary, the CSV is the evidence), the rank plot and
heatmap, `notes` — including whether the test ran on the `lognorm` snapshot
("no lognorm snapshot; markers tested on the current X" is a warning to fix).
`total_markers=0` → lower the thresholds or the domains do not differ.

**Decide.** A marker is a gene with high detection in the domain **and** low
detection in the rest (specific), not merely a large fold change. Map several
candidates per domain in space (`st_py_preprocess` spatial-plots / notebook) and
check the heatmap for cross-domain exclusivity before naming anything.

**Red lines.** Adjusted p-values are descriptive only: domains were defined from
the same data (double dipping), never report them as inferential statistics.
Tests never run on scaled values. `avg_log2FC` is this comparison's fold change,
not a portable effect size.
*Why:* Wilcoxon is the best-calibrated rank test here; detection fractions say
how specific a gene is, significance alone does not.

## 5. Domain annotation — `st_py_annotate`

**Call.** Build a panel JSON (cell type/region → gene symbols); coarse lineages
first, subtypes second. `panel_path` for scoring; `mapping_path` JSON only when
the user has ground truth. With neither, nothing is annotated (the report says
so).

**Read.** `scores[]` per cluster (per-panel scores, `genes_missing`; `best=null`
means below `min_score`), `assignments`, `label_column`, `notes` (unassigned
cluster list; absent-panel-gene warning), the dot plot and `spatial-celltypes`.

**Decide.** `min_score` starts 0; look at the winner–runner-up margin in
`scores` — a narrow win between two similar panels is not an identification,
raise `min_score` or leave the cluster unassigned. Panel symbols must match the
object's identifier type and species case. The built-in panel is a generic
example; spatial regions usually need tissue-specific markers.

**Red lines.** Zero usable panel genes = identifier mismatch, not a biology
result. Report how many clusters got each label and how many stayed unassigned —
unassigned is information, not a failure to hide. **A spatial label names the
dominant cell type or region of a spot/bin, not a single pure cell.** Never
propagate a coarse label into subtype claims.
*Why:* panel-mean-minus-rest keeps panel sizes comparable; a name must rest on
stated evidence, and a spot mixes several cells.

## 6. Spatial statistics — `st_py_spatial` (SVG + neighbourhood)

**When.** The question is spatial structure (which genes paint regions) or
spatial organisation (which domains neighbour each other). Needs the spatial
graph (§3) and, for neighbourhood, clusters.

**Call.** `svg_method=moran` (or `geary`), `svg_corr_method=fdr_bh`,
`n_perms=1000`, `cooc_interval=50`, `group_by` = the chosen cluster column.

**Read.** `svg` (`method`, `statistic`, `n_genes`, `n_significant`, `top[]` with
the statistic and `pval_adj`); `nhood_enrichment` (`top_enriched_pairs[]` with
`from`/`to`/`zscore`); `co_occurrence` (`distance_range`, `shape`); the
`svg-top`, `svg-top-gene-spatial`, `nhood-enrichment`, `nhood-cooccurrence`
figures; `notes`. A skipped section carries `{"skipped": reason}` — its absence
is never a clean result.

**Decide.** Top SVGs are candidates for domains or gradients — map them
(`svg-top-gene-spatial`) to confirm the pattern before believing biology. A
positive neighbourhood z-score marks an interface or niche (off-diagonal) or a
coherent self-adjacent region (diagonal); read it with the spatial map.
Co-occurrence above 1 at short distance = local mixing; the crossing distance is
the scale of association.

**Red lines.** Autocorrelation and enrichment depend entirely on the adjacency
graph and the cluster labels — a different `neighbor_method`/`n_neighs` or noisy
clusters changes the ranking. SVG p-values are analytic (normality assumption),
a screen, not proof of mechanism; a gene structured only by depth/batch also
scores high. Enrichment is adjacency, NOT interaction.
*Why:* Moran's I/Geary's C measure whether neighbours resemble each other;
nhood enrichment compares observed neighbour pairs to a permutation null.

## 7. Deconvolution — `st_py_deconvolve`

**When.** The question is cell-type composition of each spot/bin. **Requires a
single-cell reference** from the single-cell plugin (SingleCell agent): counts+obs
exchange or a pseudobulk table. Without a reference there is no deconvolution.

**Call.** `reference_path=<scRNA .h5ad>`, `reference_cell_type_key=cell_type`,
`deconv_method=nnls` (the only in-plugin solver).

**Read.** `deconvolution` (`method`, `reference`, `cell_types`,
`signature_genes`, `n_units`, `mean_proportions`, `note`) or `{"skipped":
reason}`; the `deconv-proportions-spatial` figure; `notes`. Per-unit proportions
are written to `obsm["deconvolution"]` and `obs["<type>_prop"]`.

**Decide.** State the reference's cell-type granularity and tissue match — a
missing type cannot be recovered, and a mismatched tissue bounds the answer.
Read the proportion maps: a type should paint a plausible region. Small
bins/cellbin are near single-cell, so deconvolution adds little there.

**Red lines.** Proportions are ESTIMATES (linear NNLS), not observed cells — never
report them as cell counts. cell2location/tangram-sc/stereopy are pip-only
project-env gaps; report them, do not install in-place. No shared genes between
object and reference = identifier mismatch, fix before fitting.
*Why:* a spot's expression is a weighted sum of cell-type profiles; NNLS solves
for non-negative weights against pseudobulk reference signatures.

## 8. Cell-cell communication — `st_py_communication`

**When.** The question is which domains might signal to each other. Needs the
spatial graph and ≥2 clusters.

**Call.** `group_by` = the chosen cluster column; `interactions_path` for a full
ligand-receptor table (TSV/CSV with source/target or ligand/receptor), or leave
it null to use the shipped curated starter (~70 pairs, offline);
`ccc_threshold=0.01`, `n_perms=1000`.

**Read.** `communication` (`method=ligrec`, `resource`,
`n_interactions_tested`, `n_significant`, `top[]` with
ligand/receptor/source_cluster/target_cluster/mean) or `{"skipped": reason}`; the
`communication-ligrec` figure; `notes`.

**Decide.** A hit is a hypothesis about an interface — confirm with the spatial
map and markers that the ligand domain and receptor domain actually neighbour
each other. The curated starter is not exhaustive; pass a fuller resource or use
liana/CellPhoneDB in a project environment for coverage.

**Red lines.** This is a SCREEN: co-expression is not proof of signalling. ligrec
runs on the lognorm snapshot (`use_raw=False`), never scaled values; results
inherit the adjacency graph and cluster labels. Gene symbols match
case-insensitively (ligrec uppercases), so a human-symbol resource works on mouse
by orthologous symbol — but verify the pairs exist in your species.
*Why:* ligrec shuffles cluster labels and asks how often mean ligand/receptor
expression exceeds the observed value, giving a permutation p-value per pair.

## 9. Cross-sample differential spatial patterns — `st_py_spatial` with `group_column`

**When.** The question is whether spatial organisation differs between conditions
(disease vs control, region A vs B). A gene can be strongly structured in one
condition and diffuse in another — invisible to ordinary differential expression.

**Call.** Set `group_column` to the condition/sample obs column; `min_units=10`.
Each condition rebuilds its OWN spatial graph and recomputes Moran's I.

**Read.** `spatial_diff` (`conditions`, `n_genes_compared`,
`top_differential_genes[]` with per-condition `moran_I` and `range`,
`condition_specific_svg`, `composition`) or `{"skipped": reason}`; the
`spatial-diff-moran` heatmap; `notes`.

**Decide.** Read the per-gene cross-condition range and each condition's SVG
list; map the pattern in each group to see what changed. Read `composition`
alongside: a shift in domain proportions can drive a shift in spatial structure.

**Red lines.** The SAMPLE is the experimental unit, not the spot — fewer than ~2
biological replicates per group → no condition-level claim; say so. This is a
descriptive screen of the difference, not a significance test of it; differences
can be driven by composition or depth shifts. Fewer than 2 usable conditions →
skipped, not a result.
*Why:* comparing the same spatial statistic across conditions isolates changed
spatial structure from changed mean expression.

## 10. The notebook path — `st_py_notebook`

**Decide.** Preset (`module` 1 QC, 2 +clustering, 3 +spatial domains/markers,
4 +SVG/neighbourhoods, 5 +deconvolution/communication, 6 full chain) or `steps`
with automatic dependency closure. Write the §2–§3 readings into `config` JSON
(`min_genes`, `min_counts`, `max_counts`, `mt_max`, `mt_pattern`, `n_hvg`,
`n_pcs`, `dims`, `resolution`, `seed`, `neighbor_method`, `n_neighs`, `platform`,
`bin_size`, `cellbin`) so the notebook and the script chain agree by
construction. `self_check=true` before delivery. `inject`/`merge` sources must
match the notebook language (python).

**Read.** The report's `steps`, `code_cells`, `self_check`; the notebook's
`phi_render` metadata (parameters with sources, fragment sha256, data
fingerprint) is the audit record.

**Red lines.** The renderer fills; you decide. Every config value carries a
source tag — a `[default]` tag on a threshold means you skipped a decision.

## 11. Red-line checklist

- Platform and analysis unit (spot / bin / cellbin) named first; Stereo-seq bin
  size chosen deliberately, not auto-accepted.
- No counts layer → re-import raw before QC/seurat_v3/deconvolution.
- No coordinates → only non-spatial analysis; SVG/nhood/communication/maps need
  the spatial adjacency graph.
- `kept_fraction` < 0.5 or zero units → stop, re-read distributions (bin too small?).
- `pct_counts_mt` absent/median 0 → wrong regex/identifier type, fix before gating.
- Markers, panel scoring, SVG and communication read `layers["lognorm"]`; never
  scaled values; "no lognorm snapshot" note → fix the chain.
- Marker p-values descriptive; fold changes not portable.
- Unassigned clusters and absent panel genes are reported, not hidden; a spatial
  label is a region/dominant type, not a pure cell.
- Spatial statistics inherit the adjacency graph — state `neighbor_method`/`n_neighs`.
- Deconvolution needs a reference; proportions are estimates, not cell counts.
- Communication/SVG are screens, not proof of signalling or mechanism.
- <2 replicates per group → no condition-level spatial claim.
- Gaps (cell2location, tangram-sc, stereopy, liana, scvi-tools…) → project
  environment, never an in-place install.

## Sources and licenses

| Source | URL | License |
| --- | --- | --- |
| OSTA — Orchestrating Spatial Transcriptomics Analysis with Bioconductor | https://bioconductor.org/books/OSTA/ | CC BY 4.0 (Crowell et al., bioRxiv 2025, doi:10.1101/2025.11.20.688607); final check pending |
| Single-cell best practices (spatial chapters) | https://www.sc-best-practices.org/ | Apache 2.0; cite Heumos, Schaar, Lance et al., Nat Rev Genet 2023, DOI 10.1038/s41576-023-00586-w |

Stereo-seq methodological points are distilled from the seven user-supplied
application papers (Dev Cell 2022, Zool Res 2024, Cell Discov 2024, Nat Genet
2024, Cell Rep 2024 ×2, Cell Metab 2024); see
`maintainers/stereo-seq-methods.md`. All prose is original, distilled from these
open sources; no upstream text, figures or data are copied. The R/Bioconductor
package names are conceptual anchors only — this plugin is the Python route
(scanpy + squidpy + spatialdata). OSTA chapter mapping:
`maintainers/osta-map.md`.