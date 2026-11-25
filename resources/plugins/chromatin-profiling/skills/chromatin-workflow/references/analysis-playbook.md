# Chromatin profiling analysis playbook (R route)

The agent's operating manual for bulk ChIP-seq and CUT&Tag/CUT&RUN. Each stage
below is a decision procedure: **Call** the tool, **Read** the named JSON fields
and figures, **Decide** the parameters from what you read, and respect the **Red
lines**. The one-line *why* gives the principle; the hard rules (pure downstream,
read-only env, lean annotation, blacklist build specificity, DiffBind-vs-csaw,
replicates as the experimental unit) live in SKILL.md `Rules`; the learner-facing
short version lives in `references/step-notes/` (one note per step). Sources and
licenses are at the end.

Quick reference: [`red-line-cheatsheet.md`](red-line-cheatsheet.md) tabulates the
hard red lines and common misreadings for §1–§9; this playbook remains the full
Call/Read/Decide procedure.

Every number here is a **starting point for ChIP/CUT&Tag data, not a default to
accept**. Read the distributions first, adapt, and say in your answer which value
you chose and why. Quote the tool's JSON numbers instead of re-deriving them. The
engine is R/Bioconductor: DiffBind 3.20 (occupancy + affinity), DESeq2/edgeR
(differential), csaw 1.44 (window-based broad marks), ChIPseeker 1.46
(annotation), ComplexHeatmap (UpSet), ggplot2 (PCA/MA). This is **pure
downstream**: it starts from called peaks, count matrices, indexed BAMs or
bigWigs — never FASTQ. Alignment and peak calling belong to the nf-core/chipseq
wrapper (ChIP) or the `nf-core/seacr/callpeak` module wrapper (CUT&Tag).

## 0. The decision loop

Chain: `chip_inspect` → `chip_qc` → `chip_consensus` → `chip_count` → then branch
on the question:

- **Which regions differ between conditions?** → `chip_diff` (DiffBind
  DESeq2/edgeR for TF/narrow marks and CUT&Tag; csaw for broad marks).
- **What genes/features do the peaks sit near?** → `chip_annotate` (needs a user
  GTF/GFF3 or TxDb).
- **Do samples separate by biology or by batch/depth?** → `chip_visualize`
  (sample PCA + peak-overlap UpSet).
- **Just explore an existing count matrix?** → `chip_visualize` with `counts`
  (skips peak loading and counting).

`chip_notebook` is the parallel teaching path; feed it the *same* threshold
readings (flags or `config` JSON). Each script tool re-runs its chain from the raw
input you give it, so pass the same input and the same thresholds every time to
keep the two paths and successive calls consistent. Skip a step only when
`chip_inspect`'s `candidate_steps` says its result already exists.

Loop at every stage: call with loose/starting values → read the JSON → decide →
re-call with the decided values → verify the counts, `warnings`, `notes`. A
default that reached the output untouched is a decision you did not make.

Call `chip_env_probe` before promising differential binding, annotation, csaw or
any optional package: the plugin environment is read-only; a gap is reported as a
fact with a project-environment next step (`env_request`), never an in-place
install. csaw on darwin-arm64, IDR, profileplyr, HOMER/MEME motif enrichment and
org.db gene symbols are gaps, not capabilities — say so, do not impersonate them.

## 1. Inspect the input — `chip_inspect`

**Call.** `chip_inspect(input, output=<project>/loaded.rds)`. With a DiffBind
sample sheet and `output`, it builds and saves the DBA object. Pure downstream —
it never expects FASTQ.

**Read.** `kind` (peaks-dir / narrowpeak / broadpeak / bed / sample-sheet /
count-table / bam-dir / bam-file / bigwig / bigwig-dir / saved-object /
annotation); `samples`; `shape` (per-sample `n_peaks`/`median_width`, or
`n_peaks`/`n_samples`/`library_sizes` for a count table, or sheet `columns`/
`conditions`/`has_bam`/`has_control`, or BAM `indexed`); `candidate_steps`;
`output` (the saved DBA path); `notes`.

**Decide.**
- **Input kind first.** A `peaks-dir`/`narrowpeak`/`broadpeak`/`bed` input can run
  qc/consensus/count/annotate. A `sample-sheet` adds `diff` (it carries Condition
  and BAM paths). A `count-table` skips straight to normalize/pca/diff. A
  `bam-dir`/`bam-file` supports qc/count/csaw. A `bigwig` is coverage-only
  (annotate/coverage-profile). A `saved-object` DBA resumes at diff/annotate.
- Whether the sheet has `bamReads` (counting and FRiP need BAM paths) and whether
  BAMs are `indexed` (a missing `.bai` blocks counting/FRiP — `samtools index` is
  upstream). Name these gaps rather than promising a step that cannot run.
- The genome build for blacklist filtering (§2/§3) and the annotation source
  (§6) — neither is bundled for the target genome.

**Red lines.** No FASTQ, no alignment, no peak calling. A count-table input has
no peaks/BAM, so FRiP and consensus are unavailable — say which input you
analysed and what it cannot support.
*Why:* the plugin is the downstream half of the workflow; upstream belongs to the
wrapper or the user's SEACR run.

## 2. Quality control — `chip_qc`

**Call.** `chip_qc(input=<peaks-dir or sample-sheet>, bam_dir=… | bam=[…],
genome=hg38|mm10 | blacklist=<BED>, figures_dir=…)`. FRiP needs indexed BAM +
matching peaks; reproducibility needs ≥2 peak sets. Pass `paired_end` correctly
(it changes how fragments are counted).

**Read.** `peaks[]` (`{sample, n_peaks, median_width}`); `frip[]` (`{sample,
total_reads, reads_in_peaks, frip, note?}`); `reproducibility`
(`overlap_fraction` matrix, `replicate_pairs[]` `{condition, a, b,
overlap_fraction}`, `count_correlation` when a count matrix exists);
`blacklist_removed`; `figures`; `warnings`.

**Decide.**
- **FRiP is the headline signal-to-noise metric.** A low FRiP means the antibody
  did not enrich — most reads are background. ENCODE guidance: TF ChIP ideally
  >1%, broad marks lower, CUT&Tag often much higher (minimal background). Read
  the per-sample `frip` distribution; a sample far below its peers is a failed
  enrichment, not a biology difference. `frip: NA` with a "no matching peak set"
  note means BAM/peak names did not align — fix the pairing, do not read it as 0.
- **Reproducibility is a proxy, not IDR.** `overlap_fraction` and
  `replicate_pairs` show whether same-condition replicates call the same peaks.
  The formal test is IDR (irreproducible discovery rate), which needs the `idr`
  pip package this environment does not ship — a reported gap. Say the overlap
  fraction is a transparent proxy, never an IDR result.
- `median_width` distinguishes narrow (TF/H3K4me3, ~200–500 bp) from broad
  (H3K27me3, kilobases) marks — this drives the DiffBind-vs-csaw choice in §5.

**Red lines.** FRiP needs indexed BAM and matching peaks; without them it is
`NA`, not a clean result. Overlap fraction is not IDR. QC runs on raw peaks/BAM,
never on normalized counts.
*Why:* FRiP measures enrichment independent of any peak caller's own report;
replicate overlap is descriptive evidence, not a calibrated reproducibility test.

## 3. Consensus peakset — `chip_consensus`

**Call.** `chip_consensus(input=<peaks-dir or sample-sheet>,
genome=hg38|mm10 | blacklist=<BED>, min_occupancy=2, figures_dir=…)`. Blacklist
filtering runs first (when a genome/blacklist is given), then the consensus is
built over the cleaned peaks.

**Read.** `n_samples`; `n_consensus`; `min_occupancy`; `occupancy_per_sample`
(how many consensus peaks each sample carries); `blacklist_removed`; `figures`;
`steps`.

**Decide.**
- **`min_occupancy` is a region-universe decision.** A site must appear in at
  least this many samples to enter the consensus — the coordinate set every
  sample is then measured over. With few replicates a high occupancy throws away
  real but variable sites; with many samples a low occupancy admits noise. Read
  `n_consensus` and `occupancy_per_sample` and say what you chose and why. The
  effective value is capped at the number of samples.
- **Blacklist before consensus.** Artifactual high-signal regions (satellite
  repeats, gaps, anomalous pileups) attract reads regardless of biology; leaving
  them in creates spurious consensus peaks. `blacklist_removed` records how many
  were dropped. Pass `genome` (hg38/mm10 inlined ENCODE v2) or `blacklist` (other
  builds); with neither, filtering is skipped and you must say so.

**Red lines.** A blacklist is genome-build specific — hg38 regions are wrong on
mm10; passing the wrong build silently removes the wrong loci. No blacklist
filtering without a genome/blacklist (the step skips and notes it).
*Why:* the consensus defines the multiple-testing universe; artifact regions
inflate it with technical signal.

## 4. Read counting — `chip_count`

**Call.** `chip_count(input=<sample-sheet or peaks-dir>, bam=[…] | bam_dir=…,
genome=… | blacklist=…, normalize_method=libsize|rpkm|none, paired_end=…,
output=<project>/counted.rds)`. Builds consensus (if not given a saved one),
counts reads per peak per sample, and normalizes for exploration.

**Read.** `n_peaks`; `n_samples`; `samples`; `library_sizes` (per-sample total
counts); `normalize_method`; `figures`; `steps`. The count matrix and normalized
matrix are written beside `output` (`counts.csv`, `normalized_counts.csv`).

**Decide.**
- **Library sizes flag depth outliers.** A sample with a vastly smaller library
  is under-sequenced; a vastly larger one may be over-amplified. Read
  `library_sizes` before trusting any between-sample comparison.
- `normalize_method`: `libsize` (CPM, default) for occupancy viz/PCA; `rpkm`
  additionally divides by peak width (matters when peaks differ greatly in
  length); `none` keeps raw counts. **This normalized matrix is for exploration
  only** — the differential engines normalize internally (§5).
- `paired_end` must match the library: single-end counts reads, paired-end counts
  fragments. A wrong setting miscounts every peak.

**Red lines.** Counting needs indexed BAM and a consensus peakset. Normalized
counts never feed the differential test — DESeq2/edgeR/csaw do their own
normalization on raw counts.
*Why:* external CPM/RPKM breaks the mean-variance relationship the
negative-binomial tests assume.

## 5. Differential binding — `chip_diff`

**Call.** `chip_diff(input=<sample-sheet>, method=deseq2|edger|csaw,
contrast=[Condition,treated,control] | (auto all pairs), min_members=2,
fdr=0.05, lfc=0, top_n=20, …)`. For `csaw`: `csaw_width=150`,
`csaw_spacing=50`, `csaw_filter=global`, `fragment_size=200`, `paired_end`.

**Read.** `method`; `fdr_threshold`; `lfc_threshold`; `n_significant`;
`n_tested` (DiffBind) or `n_windows`/`n_regions` (csaw); `table` (top regions,
FDR-sorted); `note`; `error` (a csaw platform gap surfaces here); `figures`
(MA plot); `steps`.

**Decide.**
- **DiffBind vs csaw is a biological decision.** TF and narrow marks
  (H3K4me3/H3K27ac) and CUT&Tag have sharp peaks → DiffBind (DESeq2 default,
  edgeR alternative). Broad marks (H3K27me3/H3K36me3/H3K9me3) spread over
  kilobases → csaw's sliding windows catch diffuse changes a peak-based test
  misses. Read `median_width` from §2 to confirm the mark type; say which engine
  you chose and why.
- **csaw on Apple Silicon is a gap.** bioconda has no osx-arm64 csaw build; on
  darwin-arm64 the step returns an `error` naming the gap. Use the DiffBind route
  for broad marks there, or run on linux-x64/darwin-x64, or request a project
  environment. Do not silently fall back and call it csaw.
- `method`: DESeq2 is more conservative, edgeR often more sensitive with few
  replicates — both run on the same DiffBind consensus matrix, so switching
  method does not change the peaks tested. The choice is the analyst's.
- `contrast`: an explicit `[Condition, treated, control]` tests one pair; with
  none, every condition pair with ≥`min_members` replicates is tested. csaw with
  >2 conditions tests the last level vs the first (a `note` says so) — set
  explicit contrasts for others.
- `fdr`/`lfc`: read `n_significant` vs `n_tested`. Zero significant regions with
  a clean run is a null result (or under-powered), not "no regulation"; relax
  `fdr` only with justification, and never claim significance without replicates.

**Red lines.** The sample is the experimental unit — fewer than `min_members`
(default 2) biological replicates per condition means no condition-level claim;
peaks are not replicates. Differential binding runs on raw counts with internal
normalization, never on the §4 normalized matrix. csaw needs a sample sheet with
a Condition column and indexed BAM.
*Why:* condition-level inference needs biological replicates; the NB/QL tests
assume the original count distribution and their own normalization.

## 6. Peak annotation — `chip_annotate`

**Call.** `chip_annotate(input=<peaks-dir | sample-sheet | saved .rds>,
gtf=<GTF/GFF3> | txdb=<TxDb.sqlite>, org_db=<org.Hs.eg.db|…>,
tss_upstream=3000, tss_downstream=3000, figures_dir=…)`. A saved DBA annotates its
consensus; a GRanges/GRangesList annotates directly; otherwise peaks are loaded.

**Read.** `n_peaks`; `annotation_source` (the gtf/txdb path used); `tss_region`
(`[-upstream, downstream]`); `feature_distribution[]` (`{feature, n_peaks}` —
Promoter/Exon/Intron/Intergenic…); `gene_symbols` (true only when a user org.db
is installed); `figures` (feature-distribution bar, TSS-distance); `steps`.

**Decide.**
- **Lean annotation: the target-genome TxDb comes from the user.** Pass `gtf`
  (built via `makeTxDbFromGFF`, ~a minute for a full GTF) or `txdb` (a prebuilt
  `.sqlite` via `loadDb`, faster). The plugin bundles **no** target-genome TxDb;
  without one the step reports a gap instead of guessing. The hg19 TxDb that
  ships as a ChIPseeker conda dependency is a build-time artifact, **never** the
  annotation source for hg38/mm10 peaks.
- `org_db`: gene symbols need a user org.db (e.g. org.Hs.eg.db). Without it,
  annotation still gives gene ids, features and TSS distance — `gene_symbols:
  false` means ids, not symbols. Say so; request the org.db via a project
  environment if symbols are required.
- `tss_region`: the promoter window (default ±3000 bp). Read
  `feature_distribution` — a TF peak set enriched at promoters suggests
  transcriptional regulation; a broad mark spread over gene bodies suggests
  elongation/heterochromatin. Annotation is positional inference, not measured
  target-gene regulation.

**Red lines.** No annotation without a user GTF/GFF3 or TxDb — the step reports
the gap and must not invent a target-genome annotation. The bundled hg19 TxDb is
never the default source. Annotation is inference; motif enrichment (HOMER/MEME)
is a separate gap, not part of this R downstream stack.
*Why:* ChIPseeker annotates by genomic position relative to a transcript model;
the model must match the genome the peaks were called on.

## 7. Exploration & visualization — `chip_visualize`

**Call.** `chip_visualize(input=<count-table | peaks-dir | sample-sheet>,
counts=<count.csv> (to skip counting), normalize_method=libsize|rpkm|none,
sample_column=SampleID, figures_dir=…)`. With a count table (or `counts`) it
loads counts directly; otherwise it runs the peak→consensus→count chain first.

**Read.** `n_peaks`; `n_samples`; `samples`; `normalize_method`; `pca`
(`scores[]` `{sample, condition, PC1, PC2, …}`, `variance_explained[]`,
`n_peaks`); `figures` (sample-pca, peak-upset); `steps`.

**Decide.**
- **PCA is the first check of biology vs nuisance.** Replicates of one condition
  should sit together; a PC1 split that tracks `condition` is the signal the
  differential test will quantify. A split that tracks depth/batch/antibody lot
  instead is a confound — address it before claiming differential binding. Read
  `variance_explained` to know how much PC1/PC2 carry.
- **UpSet shows peak sharing.** The overlap companion to §2 reproducibility: a
  condition whose replicates share most peaks is reproducible; samples sharing
  almost nothing differ biologically or technically.
- `normalize_method` here is for visualization only (CPM default). It does not
  feed any test.

**Red lines.** PCA/UpSet run on normalized counts for exploration; they are not
a differential test and carry no significance. A count-table input has no
occupancy matrix, so the UpSet falls back to annotated peak sets or skips.
*Why:* sample-level PCA visualizes structure; it does not model counts or
calibrate p-values.

## 8. The notebook path — `chip_notebook`

**Call.** `chip_notebook(output=<project>/chip.ipynb, data=<input>,
module=<1-6> | steps=[…], …)`. Pass the *same* thresholds read in §2–§5 as flags
or a `config` JSON (keys mirror `assets/steps.json` parameters). Pick a preset:
1 TF/narrow DiffBind+DESeq2, 2 broad csaw, 3 CUT&Tag edgeR, 4 consensus, 5
annotation, 6 count exploration; or name `steps` (dependencies are added
automatically — a `module` preset uses its explicit list, so blacklist is already
included where needed).

**Read.** The rendered notebook: one config cell at the top holding every tunable
value with a source tag, then analysis cells referencing those variables, then a
save cell. `self_check=true` executes every cell in the output directory before
delivery and fails loudly if one breaks.

**Decide.**
- `inject` adds an external R script verbatim as a code cell (e.g. a custom
  contrast or a profileplyr extension run in a project env); `merge` splices in
  another notebook whole or a cell range (`path.ipynb[0:6]`). Sources must be R —
  one kernel per notebook.
- The renderer fills placeholders from the data and the requested steps; it does
  NOT choose steps or thresholds or interpret the request — those are yours.
- Render into the user's working directory, never copy the template dir. The
  fragment library is the single analysis source: the renderer inlines fragment
  text into cells, so the notebook is self-contained in the user dir.

**Red lines.** One reading, both paths: the script chain and the notebook must
agree by construction (same flags/config), not by luck. No install cells, no
FASTQ/alignment/peak-calling cells (that is upstream). A step that needed a
missing resource says so and continues — its absence is not a clean result.
*Why:* the notebook is the teaching twin of the script chain; divergence between
them is a bug, and the read-only env cannot install.

## 9. Environment probe — `chip_env_probe`

**Call.** `chip_env_probe(module=[1-6])`. Call it before promising differential
binding, annotation, csaw or any optional package.

**Read.** `environment` (R version); `packages_total`/`packages_installed`;
`packages[]` (`{package, role, installed, version, core, required, modules[]}`);
`modules[]` (`{module, missing_required[], missing_optional[], core_ready}`);
`core_ready_modules`; `blocked_modules`; `gaps[]` (`{capability, reason,
next_step}`); `notes`.

**Decide.**
- A `core_ready` module can run here; a `blocked` module names the missing
  package and the project-environment next step. The core chain (GenomicRanges,
  rtracklayer, Rsamtools, GenomicAlignments, SummarizedExperiment, S4Vectors,
  DiffBind, DESeq2, edgeR, ChIPseeker, GenomicFeatures, ComplexHeatmap, ggplot2)
  is shipped in the managed `plugin:chip` env.
- csaw (darwin-arm64), IDR (idr pip), profileplyr, HOMER/MEME motif enrichment,
  and org.db gene symbols are gaps: report them as facts with an `env_request`
  project-environment next step, never an in-place install. The bundled
  TxDb.Hsapiens.UCSC.hg19.knownGene is a build-time ChIPseeker dependency, not a
  target-genome annotation source.

**Red lines.** The plugin environment is read-only and prebuilt — never install
into it. A gap is a fact plus a next step, not a silent skip and not an
impersonated capability.
*Why:* managed environments are reproducible and immutable; optional heavy tools
belong in a user project environment.

## 10. Red-line checklist

Before delivering, walk the table in
[`red-line-cheatsheet.md`](red-line-cheatsheet.md) (12 rows). The short version:

1. Pure downstream — called peaks/counts/BAM/bigWig in; no FASTQ, alignment or
   peak calling (ChIP → nf-core/chipseq wrapper; CUT&Tag → nf-core/seacr/callpeak module wrapper).
2. Read-only env — missing package/reference → gap report + project-env next
   step, never an in-place install or large download.
3. Lean annotation — target-genome TxDb from the user (gtf/txdb); the bundled
   hg19 TxDb is never the default source; no org.db → gene ids, not symbols.
4. Blacklist is genome-build specific — pass genome (hg38/mm10) or blacklist;
   wrong build silently removes wrong loci; neither → skip and say so.
5. DiffBind vs csaw is a biological decision — TF/narrow/CUT&Tag → DiffBind;
   broad → csaw; csaw on darwin-arm64 is a gap, not a silent fallback.
6. FRiP is the headline signal-to-noise metric, threshold from the distribution;
   IDR is a gap — overlap fraction is a proxy, not an IDR result.
7. The sample is the experimental unit — ≥min_members biological replicates per
   condition for a condition-level claim; peaks are not replicates.
8. Normalized/exploratory counts never feed the differential test — engines
   normalize internally on raw counts.
9. Consensus/min_occupancy is a region-universe decision — blacklist first, then
   choose occupancy from replicate count and biology.
10. Annotation/motif are inference/prediction, not measured regulation; motif
    enrichment (HOMER/MEME) is a gap outside this R downstream stack.
11. DiffBind needs a sample sheet + indexed BAM; relative paths resolve against
    the sheet directory; missing .bai blocks counting/FRiP.
12. The JSON is the audit record — quote its numbers; a "skipped"/gap note did
    not run, never read its absence as a clean result.

A step carrying a "skipped" note or an `error` did not run — never read its
absence as a clean result. Quote the tool's JSON numbers; do not re-derive them.

## Sources and licenses

Knowledge is an original distillation of the following open sources; no upstream
text, figures or data are copied.

| Source | URL | License |
| --- | --- | --- |
| DiffBind vignette (Stark & Brown) | https://bioconductor.org/packages/DiffBind | Artistic-2.0 |
| csaw vignette (Lun & Smyth) | https://bioconductor.org/packages/csaw | Artistic-2.0 |
| ChIPseeker vignette (Yu et al.) | https://bioconductor.org/packages/ChIPseeker | Artistic-2.0 |
| ENCODE ChIP-seq guidelines (Landt et al.) | doi:10.1101/gr.136184.111 | Genome Res 22:1813 (2012) |
| ENCODE Blacklist (Amemiya, Kundaje, Boyle) | doi:10.1038/s41598-019-45839-z | Sci Rep 9:9354 (2019); Boyle-Lab/Blacklist ENCODE v2 lists |
| nf-core/chipseq (upstream wrapper) | https://github.com/nf-core/chipseq | MIT |

The ENCODE hg38/mm10 blacklist regions are inlined (base64-gzip) in
`assets/fragments/blacklist-filter.R` so both the script tools and the rendered
notebook filter artifact regions offline. The target-genome annotation (GTF/GFF3
or TxDb) and org.db are provided by the user via `input-path` or a project
environment, never bundled for the target genome and never installed into the
read-only plugin environment. The code is original to this plugin and runs in the
private `plugin:chip` environment.