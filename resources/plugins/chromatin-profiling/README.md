# chromatin-profiling

Chromatin profiling (bulk **ChIP-seq** and **CUT&Tag/CUT&RUN**) downstream analysis plugin. **R route
only**: Bioconductor, private locked environment `plugin:chip`, tool prefix `chip_`. The main engines are
**DiffBind** (occupancy + affinity analysis) + **DESeq2/edgeR** (differential binding) + **csaw** (broad
mark window method) + **ChIPseeker** (peak annotation). The analysis unit is the **peak/sample**, not the
cell. **Pure downstream**: starts from already-called peaks/counts/indexed BAM/bigWig and never touches FASTQ.

| Input starting point | Source | Notes |
| --- | --- | --- |
| BED/narrowPeak/broadPeak | MACS2 (ChIP) / SEACR (CUT&Tag) | Preferred; per-sample peak sets |
| DiffBind sample sheet CSV | user-curated | Drives the DBA object and contrasts (Peaks/bamReads columns) |
| peak×sample count table | featureCounts / dba.count export | Skips counting; goes straight to normalize/pca/diff |
| indexed BAM | alignment output | FRiP / counting / csaw (needs .bai) |
| bigWig | coverage track | profile / visualization |
| saved `.rds` (DBA/GRanges) | this plugin's save | Resume diff/annotate |

**Alignment / peak calling is not done here**: ChIP belongs to the nf-core/chipseq wrapper; CUT&Tag upstream
SEACR has no wrapper, which is a gap, so the plugin consumes the peaks produced by the user's SEACR run directly.

## Design principles

The same contract as the single-cell / spatial / single-cell-atac plugins:

- **The plugin only does deterministic filling/rendering**. Intent parsing, step selection and threshold
  decisions belong to the agent (the rules live in
  [`agents/ChromatinProfiling.md`](agents/ChromatinProfiling.md) and the skill's SKILL.md).
- **The fragment library is the single analysis source**. The teaching code in
  `skills/chromatin-workflow/assets/fragments/` drives both the script tools (thin adapters + a JSON report
  layer) and the notebook renderer; the two forms share one copy of the code.
- **Private locked environment; do not install packages inside it**. Missing packages or missing genome
  references are reported as facts by `chip_env_probe`; the next step is a project environment (env_request)
  or instructions to download into the project directory.
- **The notebook is real code cells + per-step Markdown teaching**, rendered into the user's working
  directory, not a copied template directory.

## Three hard constraints of the R route

- **DiffBind needs a sample sheet + indexed BAM**: `dba(sampleSheet=)` builds the DBA from the CSV;
  counting/FRiP needs `.bai` indexes (`samtools index` is upstream). Relative paths resolve against the sheet
  directory.
- **csaw is a gap on Apple Silicon**: bioconda has no osx-arm64 build, so csaw only ships for
  darwin-x64/linux-x64. On darwin-arm64, broad marks take the DiffBind substitute route or move platform /
  project environment, as reported faithfully by `chip_env_probe`.
- **Figures are PNG, not plotly**: ggplot2 (PCA/MA), ChIPseeker (feature distribution/TSS distance/profile/heatmap),
  ComplexHeatmap (UpSet), rendered with ragg.

## Analysis chain

`chip_inspect` → `chip_qc` (FRiP + replicate overlap score + blacklist filtering) → `chip_consensus`
(disjoin/reduce to build the consensus peakset; min_occupancy decides the region universe) → `chip_count`
(summarizeOverlaps counting + normalization) → `chip_diff` (DiffBind DESeq2/edgeR for TF/narrow/CUT&Tag;
csaw window method for broad marks) → `chip_annotate` (ChIPseeker genomic feature distribution/nearest gene/TSS
distance + coverage profile/heatmap) → `chip_visualize` (sample PCA + peak overlap UpSet).
`chip_notebook` is the parallel teaching path; `chip_env_probe` probes before you promise heavy-dependency capabilities.

## Genome references (lean strategy)

- **blacklist inlined**: ENCODE hg38 (636 regions) / mm10 (3435 regions) is inlined as base64-gzip in
  `assets/fragments/blacklist-filter.R` and decoded offline; `genome` selects the build, `blacklist` overrides
  with another build, and with neither the step skips and says so. Blacklist regions are genome-build
  specific — hg38 regions are wrong on mm10.
- **The target-genome TxDb comes from the user**: ChIPseeker annotation uses the user's GTF/GFF3 to build a
  TxDb via `makeTxDbFromGFF`, or the user's TxDb `.sqlite` (`loadDb`); **the plugin does not bundle the target-genome
  TxDb** and reports the gap when it is missing. The hg19 TxDb shipped with the ChIPseeker conda recipe is a
  build-time dependency only and is **never the default annotation source for the target genome**.
- **Gene symbols need the user's org.db** (annoDb); without it annotation still returns gene id/features/TSS
  distance and reports the gap.
- Never install packages or download large references inside the read-only plugin environment.

## Knowledge layer

`skills/chromatin-workflow/references/` is organised by function:

- [`analysis-playbook.md`](skills/chromatin-workflow/references/analysis-playbook.md) — the agent's
  operating manual: which tool to call at each step, which JSON fields and figures to read, how to set
  thresholds, and what the red lines are.
- [`red-line-cheatsheet.md`](skills/chromatin-workflow/references/red-line-cheatsheet.md) —
  12 hard red lines × common misreadings quick-reference table, for use at runtime.
- [`step-notes/`](skills/chromatin-workflow/references/step-notes/) — the purpose/principle/interpretation/pitfalls
  of each step, i.e. the notebook's teaching cells (19 steps).
- [`maintainers/implementation-notes.md`](skills/chromatin-workflow/references/maintainers/implementation-notes.md)
  — implementation notes and provenance (fragment library architecture, DiffBind/csaw/ChIPseeker API mapping,
  lean annotation and blacklist self-containment, JSON contract).

The knowledge is an original distillation of the [DiffBind](https://bioconductor.org/packages/DiffBind),
[csaw](https://bioconductor.org/packages/csaw), [ChIPseeker](https://bioconductor.org/packages/ChIPseeker)
vignettes (Artistic-2.0) and the ENCODE ChIP-seq guidelines (doi:10.1101/gr.136184.111); the blacklist comes
from doi:10.1038/s41598-019-45839-z. No upstream text/figures/data are copied.

## Environment

`environments/chip/environment.yml` declares the dependencies (pure conda + R source packages, no pip), and
`locks/` holds the lock files for the three platforms (darwin-arm64, darwin-x64, linux-x64). r-base is pinned
to 4.5 (DiffBind 3.20 / ChIPseeker 1.46 / csaw 1.44 / DESeq2 1.50 / edgeR 4.10, the same line as the
single-cell seurat route). csaw is a platform-conditional dependency (CondaPlatformDependency, excluding
darwin-arm64). GenomeInfoDbData goes through a GitHub mirror pin (codeload reachable). Solving goes through the
tuna mirror, then the host is rewritten back to the official one.

## Gaps (reported as facts by `chip_env_probe`, handled via a project environment)

- **csaw (darwin-arm64)**: bioconda has no osx-arm64 build; on Apple Silicon broad marks go the DiffBind
  substitute route or move to linux-x64/darwin-x64/project environment.
- **IDR**: the formal reproducibility test needs the `idr` pip package, and the locked environment does not
  support pip; the plugin provides a transparent overlap-score proxy and does not pass it off as IDR.
- **profileplyr**: r44 only, incompatible with the r-base 4.5 main line, not included.
- **motif enrichment (HOMER/MEME)**: not in the R downstream stack, a gap; motif hits are sequence
  predictions, not measured TF binding.
- **org.db gene symbols**: needs the user's org.db; when missing, gene ids are reported.
- **CUT&Tag upstream SEACR**: no wrapper; the plugin consumes the user's SEACR peaks; alignment/peak calling
  are an upstream gap.

## Relationship to other plugins

- **single-cell-atac** (sibling plugin): scATAC is another plugin (Python/snapatac2, cell-level accessibility),
  kept separate from this plugin by analysis paradigm — this plugin is bulk, peak/sample-level, R route.
- **wrappers**: raw reads to BAM/peak (alignment, MACS2) is the nf-core/chipseq wrapper's job, not a notebook cell.
- **single-cell**: if an scRNA reference is needed (peak annotation against expression), reference objects come
  via the SingleCell agent's shared exit.