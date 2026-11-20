# Read-based profiling playbook

The decision loop for every step: **call → read the JSON → decide → check red
lines**. Tools are deterministic; you (the agent) choose steps and thresholds.
Upstream classification runs through wrappers; this chain starts at reports.

## 0. Route check

- Question is "who is there / what are they doing / how do groups differ" →
  this skill. "Recover genomes / novel species / strain resolution" →
  `mag-workflow`. Both can run on the same reads; they answer different
  questions.
- 16S/amplicon data does not belong here — use the `ampliseq` workflow
  wrapper. This skill is for shotgun metagenomes.

## 1. Upstream (wrappers, not plugin scripts)

Typical chains (all via `wrapper_search` → `wrapper_inspect` →
`wrapper_plan_run` → `wrapper_run`):

- QC + host removal: `fastp` → `kneaddata` (local wrapper; host database
  required) or `bowtie2` against a host index. EasyMetagenome defaults:
  fastp qualified quality 15, length required 50; subsample to ~20M reads
  (80M lines) with a fixed seed for cross-sample comparability when depth
  varies wildly.
- k-mer route: `kraken2` (database: standard ~50GB+, confirm disk first) →
  `bracken` (read length must match the library, threshold default 10) →
  optional `krakentools` kreport2mpa / `krona`.
- marker-gene route: `metaphlan` (v4 database ~5GB; 6-8 cores / 32GB RAM
  recommended by the cambiotraining course) → optional `humann3` for
  function (HUMAnN runs MetaPhlAn internally; ~15GB databases).
- Databases are the hard constraint. Confirm paths and free disk before any
  run; ask the user before downloading. Never download silently.

## 2. Inspect — `mg_inspect`

- Call on the report directory or file. Read `files[].kind`, `samples`,
  `ready_steps`.
- Decide: one kind → merge; mixed kinds → analyze separately; `unknown` →
  open the file yourself (read tool) before guessing.
- Red line: do not assume value kind. Kraken2/Bracken are counts; MetaPhlAn
  is relative abundance (0–100). The merge result states it; carry it forward.

## 3. Merge — `mg_merge`

- Call with all per-sample reports of one kind, `kind`, and a rank.
- Rank choice: species (S) for composition/differential questions; genus (G)
  or phylum (P) when species-level noise dominates (sib-swiss course works at
  phylum for ordination and warns species plots get messy).
- Read `taxa_before`/`taxa_after`. A rank with near-zero taxa means the wrong
  rank code or a database without that resolution — stop and check.
- Filters: `min_abundance` removes rare noise; `top` is for figures, not for
  diversity (diversity needs the full community). Never pass a top-filtered
  matrix to `mg_diversity`.
- Red line: one matrix = one kind + one rank + one value kind.

## 4. Rarefaction — `mg_rarefaction`

- Call on the merged **count** matrix before comparing richness across
  samples; pass metadata to color curves by group.
- Read `library_sizes` and `not_plateaued`: samples still rising at their own
  depth have depth-limited richness — say so before quoting observed counts.
  The dashed line marks the smallest library, the fair common depth.
- Decide: all curves plateau well below the smallest library → observed
  richness is comparable; otherwise report Chao1 alongside, or compare at a
  common depth.
- Red lines: counts only — the tool rejects MetaPhlAn relative abundances and
  fractional Bracken estimates. Rarefaction discards data: it is a depth
  check, never a preprocessing step for `mg_differential`.

## 5. Diversity — `mg_diversity`

- Call with the full matrix, correct `values`, metadata + `group_col` when
  groups exist.
- Read: alpha table (Shannon comparable across samples only at similar
  depth — check `mg_rarefaction` first if depths differ), PCoA variance
  explained (first two axes < ~30% means the 2D picture is weak; say so),
  PERMANOVA pseudo-F and p (seeded; report permutations).
- Decide: group separation on PCoA + PERMANOVA p < 0.05 → composition
  differs; then run differential. No separation → do not fish for taxa.
- Red lines: Chao1/observed on relative abundances are invalid (the tool
  skips them — report the skip, do not impute). Bray-Curtis on CLR-transformed
  data is a different analysis; if the user wants Aitchison geometry, that is
  a project-environment task (microbiome/ANCOM-BC are not in this plugin).

## 6. Differential — `mg_differential`

- Call with matrix + metadata; name `case`/`control` explicitly for two-group
  contrasts; >2 groups gives Kruskal-Wallis (omnibus; post-hoc is your job).
- Read `n_significant`, top taxa, volcano. Defaults FDR 0.05, |log2FC| ≥ 1
  are starting points — justify any change.
- Decide: significant taxa drive the biological story; check they are not
  all one clade at a shallow rank (aggregate to genus/family and re-test if
  species calls are unstable).
- Red lines: Wilcoxon on dozens of samples with many zero-inflated taxa has
  low power — report the limitation. For LEfSe-style LDA scores and
  cladograms, run the `lefse` wrapper then `graphlan` wrapper on the same
  matrix (EasyMetagenome §2.3 / Snakemake Meta A track); this tool does not
  fake LDA.

## 7. Function — `mg_function`

- HUMAnN genefamilies are stratified (gene|bug). `unstratified` gives the
  community-total table for group comparisons; `stratified` keeps
  per-species contributions for attribution.
- `normalize` (relative/CPM) before cross-sample comparison; `regroup` with
  a feature→group mapping (e.g. UniRef90→KO, EC→pathway) to change the
  functional level. EasyMetagenome §2.4 regroups to KO then draws GraphLan.
- Red line: this is plain column renormalization, not HUMAnN's
  stratification-aware renorm — for exact HUMAnN semantics use the humann3
  wrapper utilities. Do not compare stratified sums to unstratified totals
  as if they were the same quantity.

## 8. Figures — `mg_plot`

- Barplot top 20 (default) at the rank you analyzed; group-mean barplot when
  metadata exists. Heatmap with `trans: log10` for counts spanning orders of
  magnitude; add `cluster: rows` (or `both`) to group co-occurring taxa and
  similar samples. Clustering reorders only — dendrograms and publication
  palettes are the visualization plugin's job.
- Red line: figures show what the matrix holds — never relabel relative
  abundance as reads.

## 9. Red-line checklist (whole chain)

- [ ] Database paths/disk confirmed before any upstream wrapper run.
- [ ] Value kind carried explicitly from merge to every consumer.
- [ ] Depth adequacy checked (rarefaction) before quoting observed richness.
- [ ] Diversity ran on the full matrix, not a top-N subset or rarefied data.
- [ ] Thresholds stated with reasons; JSON counts cited, not recomputed.
- [ ] Gaps named: LEfSe/cladogram, Krona, phyloseq networks, ANCOM-BC,
      StrainPhlAn — wrappers or project environment, never pretended.

## Sources and licenses

See [maintainers/source-map.md](maintainers/source-map.md). Course ideas are
distilled with attribution; no GPL code was copied into this plugin.