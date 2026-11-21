# Implementation notes — chromatin-profiling R route

What each tool computes and how the layers divide the work. The chain is
deterministic: every step reads the shared bookkeeping, updates it, and reports the
parameters it used. The analysis lives in `assets/fragments/` (per-step teaching
notes in `references/step-notes/`) — the tools are thin adapters over the same
fragments the notebook renders, so **adapt the fragment, not the script**. The
engine is R/Bioconductor: DiffBind 3.20 (occupancy + affinity), DESeq2/edgeR
(differential), csaw 1.44 (window-based broad marks), ChIPseeker 1.46 (annotation),
GenomicFeatures (TxDb), GenomicAlignments/Rsamtools (counting), ComplexHeatmap
(UpSet), ggplot2 (PCA/MA). This is **pure downstream**: called peaks, count
matrices, indexed BAMs or bigWigs in — never FASTQ.

## Architecture

```
assets/steps.json            registry: parameters (name/type/default/comment),
                             steps (id/title/requires), presets, kernel, footer
assets/fragments/<step>.R    self-contained teaching code, sys.source'd into
                             globalenv; the SINGLE analysis source for both paths
references/step-notes/<step>.md  the notebook's per-step Markdown teaching cell
scripts/_common.R            flag parsing (parse_flags/flag_*), emit/fail/info/chr,
                             require_bioc, read_object, import_peaks,
                             read_count_table, save_plot, script_dir
scripts/_config.R            builds EVERY parameter global from steps.json defaults
                             + flag overrides; adds structural data_path/output/
                             figures_dir; auto-detects a sample sheet from data_path
scripts/_fragments.R         run_fragments(steps) sys.sources fragments in order
                             into globalenv; nothing analytical lives here
scripts/<tool>.R             thin adapter: parse flags → require packages →
                             sink(stderr()) → run_fragments → sink() → report;
                             stdout = exactly one JSON
scripts/reports/<name>.R     serializes fragment bookkeeping into the JSON contract
scripts/env_probe.R          package/module probe (read-only env facts)
schemas/env-probe.json       the one declared tool output schema (env_probe only)
scripts/build_notebook.py    renderer (RENDERER_VERSION=3, language-neutral): reads
                             fragments + step-notes, emits config cell + code/MD cells
```

**Why `_config.R`.** A fragment references ~26 parameter globals (`min_occupancy`,
`fdr_threshold`, `csaw_width`, `tss_upstream`, …). Building them from the registry
guarantees every global exists with its declared default, then overridden by flags.
This keeps the adapter namespace a superset of the notebook config namespace — the
two paths agree by construction, and a fragment never hits an "object not found"
because a flag was omitted. **CLI flag names diverge from steps.json param keys**:
`_config.R` reads `peaks/counts/bam/sample_sheet/blacklist/fdr/lfc/gtf/txdb`, while
steps.json params are `peak_paths/count_path/bam_paths/sample_sheet_path/
blacklist_path/fdr_threshold/lfc_threshold/gtf_path/txdb_path`. SKILL analysis-tool
args use the CLI names; the notebook tool args use the param keys (the renderer
auto-generates `--<param>` flags). `parse_flags` normalizes dashes to underscores.

**stdout discipline.** Each adapter wraps `run_fragments` in `sink(stderr())` so
fragment chatter (`cat`, package messages) goes to stderr; the report's `emit`
(`jsonlite::toJSON(auto_unbox=TRUE, null="null", digits=NA, pretty=TRUE)`) writes
exactly one JSON object to stdout. A tool's stdout is always parseable JSON.

**Bookkeeping.** `setup.R` initializes the shared globals; fragments leave
`figure_paths`, `analysis_steps`, `analysis_notes`, `peak_list`, `consensus_peaks`,
`occupancy_matrix`, `count_matrix`, `normalized_counts`, `dba_object`, `csaw_result`,
`diff_results`, `qc_metrics`, `blacklist_removed`, `annotation_result`, `txdb`,
`sample_sheet`, `pca_result` in globalenv. Reports read these; they never recompute
analysis. `show_plot` prints inline only under IRkernel, otherwise writes a ragg PNG
into figures_dir and records the path.

## R-specific constraints

- **Top-level `if/else` parse gotcha.** `x <- if (cond) a` newline `else b` FAILS at
  top level under `sys.source` ("unexpected 'else'"). Use braces `{ }` or put `else`
  on the same line. Inside braces/function bodies/list() it is fine. The annotate,
  coverage-profile and reports/annotate blocks use brace form for this reason.
- **DiffBind needs a sample sheet.** `dba(sampleSheet=)` builds the DBA from a CSV
  with SampleID/Condition/Replicate/Peaks/PeakCaller (required); bamReads is needed
  for counting/FRiP, ControlID/bamControl optional. Relative Peaks/bamReads paths
  resolve against the sheet directory (load-peaks/load-bam do this).
- **csaw is platform-gated.** bioconda has no osx-arm64 build; `CondaPlatformDependency`
  excludes csaw from darwin-arm64 in the env lock. `csaw.R` checks `have("csaw")` and
  returns an `error` naming the gap rather than aborting — the report surfaces it.
- **Lean annotation.** No target-genome TxDb is bundled. `annotate.R`/`coverage-profile.R`
  load a user `txdb_path` (AnnotationDbi::loadDb) or build from `gtf_path`
  (txdbmaker::makeTxDbFromGFF — the builder split out of GenomicFeatures in Bioc >= 3.22);
  without either they stop with a gap. The hg19 TxDb that ships as
  a ChIPseeker conda dependency is build-time only, never the annotation source. Gene
  symbols need a user org.db (annoDb); without it, gene ids are reported.
- **Blacklist self-containment.** ENCODE v2 hg38 (636 regions) / mm10 (3435 regions)
  are inlined base64-gzip INSIDE `blacklist-filter.R` (not separate BED assets),
  because the rendered notebook runs in the user dir with no plugin asset tree.
  Decode: `memDecompress(jsonlite::base64_dec(b64), type="gzip", asChar=TRUE)` →
  `cat(bed, file=tempfile)`. `blacklist_path` overrides; `genome` selects the build;
  neither → skip. Same payloads as the scatac plugin's blacklist.py.
- **Figures are PNG**, not plotly: ggplot2 (PCA/MA), ChIPseeker (feature bar, TSS
  distance, profile, heatmap), ComplexHeatmap (UpSet). `save_plot`/`show_plot` use
  ragg when available.

## Input loading

`project_inspect.R`'s `detect_kind` classifies the input; the load fragments
dispatch. Pure downstream — it never expects FASTQ:

| kind | loader | notes |
| --- | --- | --- |
| peaks-dir / narrowpeak / broadpeak / bed | `load-peaks` (rtracklayer::import → GRangesList) | per-sample n_peaks/median_width; names = file basenames |
| sample-sheet | `load-sample-sheet` (DiffBind::dba) | builds DBA; Peaks/bamReads drive load-peaks/load-bam |
| count-table | `load-counts` (read.delim → numeric matrix) | no peaks/BAM → FRiP/consensus unavailable |
| bam-dir / bam-file | `load-bam` (Rsamtools::BamFileList) | needs .bai; missing index warned |
| bigwig / bigwig-dir | coverage-only | annotate/coverage-profile |
| saved-object (.rds) | `read_object` | DBA → diff/annotate; GRanges → annotate/count |
| annotation (GTF/GFF3) | TxDb source | annotate only |

## inspect

`project_inspect.R` runs `detect_kind` and reports `kind`, `samples`, `shape`,
`candidate_steps`, `notes`; with a sample sheet + `output` it builds and saves the
DBA (`DiffBind::dba(sampleSheet=)` → `write_object`). It requires GenomicRanges +
rtracklayer; DiffBind only when building a DBA. No fragments run here — it is a
read-only classifier (plus optional DBA construction).

## qc

`qc_filter.R` runs `setup → [load-sample-sheet] → load-peaks → load-bam →
blacklist-filter → qc-frip → qc-reproducibility [+ save]`. `qc-frip.R` uses
Rsamtools::countBam (mapped reads) + GenomicAlignments::summarizeOverlaps (Union) per
sample, pairing BAM to peaks by name then position; `qc-reproducibility.R` builds the
overlap_fraction matrix (overlapsAny), same-condition replicate_pairs from the sheet,
and a Spearman count_correlation when a count matrix exists. Report emits `peaks[]`,
`frip[]`, `reproducibility`, `blacklist_removed`, `figures`, `warnings` (IDR-is-a-gap
note always present).

## consensus

`consensus_peaks.R` runs `setup → [load-sample-sheet] → load-peaks → blacklist-filter
→ consensus [+ save]`. `consensus.R` disjoins all peaks, keeps intervals present in
≥min_occupancy samples (capped at n_samples), reduces to consensus_peaks, and
recomputes occupancy_matrix on the merged set. Report emits `n_samples`,
`n_consensus`, `min_occupancy`, `occupancy_per_sample`, `blacklist_removed`.

## count

`count_reads.R` runs `setup → [load-sample-sheet] → load-peaks → load-bam →
blacklist-filter → consensus → count → normalize [+ save]`. `count.R` uses
summarizeOverlaps(consensus, bam_files, Union, inter.feature=FALSE) → count_matrix,
and dba.count(DBA_SCORE_READS) when a DBA exists. `normalize.R` computes CPM/RPKM/none
→ normalized_counts and dba.normalize(normalize = DBA_NORM_LIB). Note `dba.normalize`'s
`method=` selects the differential engine, not the normalization scheme — the scheme
token (DBA_NORM_LIB) goes to `normalize=`. Report emits `n_peaks`,
`n_samples`, `samples`, `library_sizes`, `normalize_method`.

## diff

`diff_binding.R` validates `method` in deseq2/edger/csaw. csaw route: `setup →
load-sample-sheet → load-bam → csaw`; else `setup → load-sample-sheet → load-peaks →
load-bam → blacklist-filter → consensus → count → diff-deseq2|diff-edger [+ save]`.
`diff-deseq2.R`/`diff-edger.R`: dba.contrast (explicit contrast or DBA_CONDITION,
minMembers) → dba.analyze(DBA_DESEQ2|DBA_EDGER) → dba.report(method=DBA_DESEQ2|DBA_EDGER,
th=fdr, fold=lfc, bUsePval=FALSE) → diff_results (method/n_significant/n_tested/table) +
dba.plotMA(method=…). `dba.report` defaults to DESeq2, so `method=` MUST match the engine
`dba.analyze` just ran — otherwise an edgeR run fails with "DESeq2 analysis has not been
run for this contrast".
`csaw.R`: windowCounts → filterWindows → normFactors (10kb bins) → asDGEList →
estimateDisp/glmQLFit/glmQLFTest → mergeWindows/combineTests → csaw_result
(n_windows/n_regions/n_significant/table); a platform gap returns an error. Report
picks `csaw_result %||% diff_results`, emits `method`, `n_significant`, `n_tested`/
`n_windows`/`n_regions`, `table`, `note`, `error`.

## annotate

`annotate_peaks.R` runs `setup` FIRST, then routes a saved .rds (DBA →
dba.peakset(bRetrieve=TRUE) → consensus_peaks; GRanges → consensus_peaks;
GRangesList → peak_list) or `[load-sample-sheet] → load-peaks`, then `annotate →
coverage-profile [+ save]`. `annotate.R` loads/builds the user TxDb, then annotatePeak
per sample (tssRegion, level="gene", annoDb when org.db installed) → annotation_result
(a csAnno for a consensus GRanges, or a per-sample list via lapply), plotAnnoBar +
plotDistToTSS. `coverage-profile.R` plotAvgProf2(peaks) + peakHeatmap(peaks) (the tag
matrix is built internally). Report handles csAnno vs list (do.call rbind), emits
`n_peaks`, `annotation_source`, `tss_region`, `feature_distribution`, `gene_symbols`.

## visualize

`visualize.R` branches: a count table (or `counts`) → `load-counts`; else
`[load-sample-sheet] → load-peaks → load-bam → blacklist-filter → consensus → count`.
Then `normalize → pca → upset [+ save]`. `pca.R` log2(normalized+1) over variable
peaks → prcomp → pca_result (scores/variance_explained/n_peaks) + ggplot colored by
Condition. `upset.R` ComplexHeatmap::make_comb_mat(occupancy_matrix) → UpSet, else
ChIPseeker::upsetplot on a single csAnno when ggupset is installed, else skip. Report emits `n_peaks`,
`n_samples`, `samples`, `normalize_method`, `pca`, `figures`.

## env_probe

`env_probe.R` probes a 20-package catalogue across 6 modules via requireNamespace,
reports present/absent as facts with role/core/required/modules. Core required:
GenomicRanges/rtracklayer/Rsamtools/GenomicAlignments/SummarizedExperiment/S4Vectors/
DiffBind/DESeq2/edgeR/ChIPseeker/GenomicFeatures/ComplexHeatmap/ggplot2; ragg optional.
Gaps: csaw (darwin-arm64), idr, profileplyr, org.db (when neither Hs/Mm present),
HOMER/MEME motif (static), upstream peak calling (static). The bundled
TxDb.Hsapiens.UCSC.hg19.knownGene is listed optional with modules=[] and role
"build-time ChIPseeker conda dependency; NOT a target-genome annotation source". The
one declared tool output schema is `schemas/env-probe.json` (draft 2020-12,
additionalProperties:false). Absent packages are project-environment next steps
(`env_request`), never in-place installs (the plugin env is read-only, conda-locked).

## Environment lock ordering

The `environments/chip/locks/*.txt` files are NOT in Phi's native solve order — they are
**stably re-partitioned into three groups, in this order**: (1) all **conda-forge** URLs,
(2) **bioconda** regular packages, (3) **bioconda data packages** (go.db, org.*.db, txdb.*,
kegg.db, reactome.db, BSgenome.*). Relative order within each group is preserved, and a
`# order: …` header comment records the rule.

Why: Phi's `lock-env.ts` emits a URL-alphabetical `@EXPLICIT` lock and `micromamba create
-f` links **in file order**. Two independent post-link hazards follow from that order:

1. `conda.anaconda.org/bioconda/…` sorts before `conda-forge/…`, so `yq` and `r-base`
   (conda-forge) link *after* the Bioconductor data packages. Their post-link
   `installBiocDataPackage.sh` needs `yq` (JSON parsing) and `Rscript` / `R CMD INSTALL` →
   "yq: command not found". micromamba 2.1.1 only warns (sqlite never materialised → the
   package silently breaks later); micromamba 2.9.0 adds a
   `libmamba Executing post-link script …` marker that Phi's `link-script-failures.ts` turns
   into a hard error.
2. Within bioconda, `bioconda/noarch/…` sorts before `bioconda/<platform>/…`, so the noarch
   data packages (go.db, txdb.*) link *before* the platform-specific Bioconductor runtime
   packages they depend on (biobase, S4Vectors). The data-package post-link's `R CMD INSTALL`
   then fails with `package 'Biobase' required by 'AnnotationDbi' could not be found` /
   `package 'S4Vectors' required by 'GenomicFeatures' could not be found`, and the sqlite is
   never installed.

Putting all data packages last (after every other bioconda package) guarantees both `yq`/
`r-base` and every R dependency are present when their post-link runs. bioconda depends on
conda-forge (never the reverse), so conda-forge-first is safe. The reorder is idempotent
(restores the order note) and **must be re-applied after every re-solve**, because Phi's tool
always re-emits URL order. The env spec also carries `bioconductor-txdbmaker` (the GTF→TxDb
builder split out of GenomicFeatures in Bioc ≥ 3.22) and `r-ggupset` (ChIPseeker::upsetplot
Suggests; the single-csAnno UpSet fallback); both are required by the fragments.

## Bioconductor API mapping (quick reference)

| step | call | result |
| --- | --- | --- |
| load-peaks | `rtracklayer::import` | GRangesList peak_list |
| load-sample-sheet | `DiffBind::dba(sampleSheet=)` | dba_object |
| load-bam | `Rsamtools::BamFileList` | bam_files |
| load-counts | `read.delim` | count_matrix |
| qc-frip | `Rsamtools::countBam`, `GenomicAlignments::summarizeOverlaps` | qc_metrics$frip |
| qc-reproducibility | `GenomicRanges::overlapsAny`, `stats::cor` | qc_metrics$reproducibility |
| blacklist-filter | `memDecompress` + `rtracklayer::import` + subset | blacklist_removed |
| consensus | `GenomicRanges::disjoin/reduce/overlapsAny` | consensus_peaks, occupancy_matrix |
| count | `summarizeOverlaps`, `DiffBind::dba.count` | count_matrix |
| normalize | CPM/RPKM, `DiffBind::dba.normalize` | normalized_counts |
| pca | `stats::prcomp`, ggplot2 | pca_result |
| diff-deseq2/edger | `dba.contrast/analyze/report/plotMA` | diff_results |
| csaw | `windowCounts/filterWindows/normFactors`, edgeR QL, `mergeWindows/combineTests` | csaw_result |
| annotate | `AnnotationDbi::loadDb` / `txdbmaker::makeTxDbFromGFF`, `ChIPseeker::annotatePeak` (per-sample `lapply` → csAnno list) | annotation_result, anno_df, txdb |
| coverage-profile | `ChIPseeker::plotAvgProf2/peakHeatmap` (take peaks directly; tag matrix built internally) | figures |
| upset | `ComplexHeatmap::make_comb_mat/UpSet` | figures |
| save | `saveRDS`, `rtracklayer::export`, `write.csv` | .rds + BED/CSV tables |

## Preset mapping

| Module | Preset | Chain coverage |
| --- | --- | --- |
| 1 TF/narrow ChIP | `module: 1` | setup → load-sample-sheet → load-peaks → load-bam → blacklist-filter → consensus → count → normalize → diff-deseq2 → save (15 steps) |
| 2 Broad mark | `module: 2` | setup → load-sample-sheet → load-bam → csaw → save (10 steps) |
| 3 CUT&Tag | `module: 3` | DiffBind chain with diff-edger (15 steps) |
| 4 Consensus | `module: 4` | setup → load-peaks → blacklist-filter → consensus → save (6 steps) |
| 5 Annotation | `module: 5` | setup → load-peaks → annotate → coverage-profile → save (5 steps) |
| 6 Count exploration | `module: 6` | setup → load-counts → normalize → pca → upset (5 steps) |

A `module` preset uses its explicit step list with NO dependency closure, so
blacklist is listed explicitly where needed; `steps` closes over `requires`. Steps
needing a missing package/reference report the gap and continue (the notebook footer
lists what did not run). FASTQ→BAM→peak (nf-core/chipseq for ChIP; user SEACR for
CUT&Tag, a wrapper gap) is upstream, not notebook cells. scATAC is the sibling
`single-cell-atac` plugin (Python/snapatac2).

## JSON contract

Every tool emits one object with `input`/`output` plus its section blocks. Common
conventions: `figures[]` lists PNG paths actually written; `steps[]` lists the
fragments that ran; `warnings`/`notes`/`error` carry gaps and skipped-step reasons —
a skipped step appends a note or an `error`, never an empty success. Only `env_probe`
declares a tool-level output schema (`schemas/env-probe.json`); the other tools'
contracts are documented here and in the playbook. Adapter and notebook share the
fragment library, so the same reading fed to both produces the same numbers by
construction. The renderer (`build_notebook.py`, RENDERER_VERSION=3, language-neutral)
reads each `.R` fragment plus `references/step-notes/<step>.md` (a missing note →
empty Markdown cell, code still runs), builds one config cell from steps.json
parameters + flags/config JSON (precedence flag > config > default, `<-` assignment,
source tags), and ends with a save cell. `--inject`/`--merge` sources must be R (one
kernel per notebook).