---
name: chromatin-workflow
description: >-
  Run a bulk ChIP-seq / CUT&Tag / CUT&RUN downstream workflow on called peaks
  (BED/narrowPeak/broadPeak), a DiffBind sample sheet, a peak x sample count table,
  indexed BAMs, a bigWig or a saved .rds: inspect, QC (FRiP + replicate
  peak-overlap), blacklist filtering, consensus peakset, counting + normalization,
  PCA + UpSet, differential binding (DiffBind/DESeq2/edgeR for TF/narrow marks +
  CUT&Tag; csaw window-based for broad marks), and ChIPseeker annotation / TSS
  profiles against a user GTF/GFF3 or TxDb. Use for ChIP-seq, CUT&Tag, 染色质
  profiling, 组蛋白修饰, 转录因子结合, DiffBind, csaw, ChIPseeker, FRiP, 差异结合,
  broad mark, H3K27me3, H3K4me3, peak 注释, .narrowPeak, .broadPeak. R/Bioconductor
  route (DiffBind 3.20 + ChIPseeker 1.46 + csaw 1.44). Pure downstream:
  alignment/peak calling are the nf-core/chipseq wrapper (ChIP) or
  nf-core/seacr/callpeak module wrapper (CUT&Tag). csaw on darwin-arm64, IDR, profileplyr, HOMER/MEME motif
  and org.db gene symbols are reported gaps. Not for single-cell ATAC (the
  single-cell-atac sibling).
phi:
  environment: plugin:chip
  attachTo:
    - ChromatinProfiling
  scripts:
    - name: inspect
      description: >-
        Read a chromatin-profiling input and report its kind (peaks file/dir,
        DiffBind sample sheet, count table, BAM file/dir, bigWig, saved .rds or
        annotation GTF), shape, per-sample peak counts and median width, sample
        sheet columns and conditions, BAM index status, and which follow-up steps
        the input is ready for. Call this first on an unfamiliar input. With
        output and a sample sheet, builds and saves a DiffBind DBA object so later
        steps can chain.
      run: [Rscript, ./scripts/project_inspect.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: Called peaks (BED/narrowPeak/broadPeak file or directory), a DiffBind sample sheet CSV, a peak x sample count table, BAM file/dir, bigWig, a saved .rds, or an annotation GTF/GFF3.
          output:
            type: string
            format: project-path
            description: Where to save a DiffBind DBA object as .rds (sample-sheet input only), inside the project.
      approval: read
      output: ./schemas/step-result.json
      timeoutSeconds: 600
    - name: qc
      description: >-
        Compute per-sample ChIP/CUT&Tag quality metrics from called peaks plus
        indexed BAMs: FRiP (fraction of reads in peaks, the headline enrichment
        metric) and replicate reproducibility (pairwise peak-overlap fractions and
        same-condition replicate pairs, with a count correlation when counts
        exist). Optionally removes ENCODE blacklist regions first. Reports the
        metrics behind a sensible enrichment judgment; IDR is the formal
        reproducibility test and is a reported gap (overlap fraction is a proxy).
        Start from peaks + BAMs or a sample sheet with bamReads.
      run: [Rscript, ./scripts/qc_filter.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A peaks directory/file, or a DiffBind sample sheet CSV with Peaks and bamReads columns.
          output:
            type: string
            format: project-path
            description: Where to save the QC'd object/peaks as .rds.
          figures_dir:
            type: string
            format: project-path
            description: Directory for QC figures; omit to skip figures.
          peaks:
            type: array
            items:
              type: string
              format: input-path
            description: Explicit peak files; empty derives them from the sample sheet Peaks column or input.
          bam:
            type: array
            items:
              type: string
              format: input-path
            description: Explicit indexed BAM files; empty derives them from the sample sheet bamReads column, bam_dir or input.
          bam_dir:
            type: string
            format: input-path
            description: Directory of .bam files (with .bai indexes).
          sample_sheet:
            type: string
            format: input-path
            description: DiffBind-format CSV; null uses input when it is a sample sheet.
          paired_end:
            type: boolean
            description: Sequencing is paired-end (default false); affects FRiP counting.
          blacklist:
            type: string
            format: input-path
            description: User ENCODE blacklist BED; overrides the inlined hg38/mm10 lists.
          genome:
            type: string
            description: Genome build (hg38 or mm10) selecting the inlined ENCODE v2 blacklist; other builds need blacklist.
          sample_column:
            type: string
            description: Sample sheet column naming each sample (default SampleID).
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 1800
    - name: consensus
      description: >-
        Build a reproducible consensus peakset from called peaks: optionally remove
        blacklist regions, disjoin every sample's peaks into non-overlapping
        intervals, keep intervals present in at least min_occupancy samples, merge
        adjacent kept intervals, and visualise sample overlap with an UpSet plot.
        This is the region universe counting and differential binding measure over.
        Reports the consensus size and per-sample occupancy. Start from peaks or a
        sample sheet.
      run: [Rscript, ./scripts/consensus_peaks.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A peaks directory/file, or a DiffBind sample sheet CSV with a Peaks column.
          output:
            type: string
            format: project-path
            description: Where to save the consensus GRanges (and tables) as .rds.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the overlap UpSet figure.
          peaks:
            type: array
            items:
              type: string
              format: input-path
            description: Explicit peak files; empty derives them from the sample sheet or input.
          sample_sheet:
            type: string
            format: input-path
            description: DiffBind-format CSV; null uses input when it is a sample sheet.
          blacklist:
            type: string
            format: input-path
            description: User ENCODE blacklist BED; overrides the inlined hg38/mm10 lists.
          genome:
            type: string
            description: Genome build (hg38 or mm10) selecting the inlined ENCODE v2 blacklist.
          min_occupancy:
            type: integer
            minimum: 1
            description: A consensus peak must appear in at least this many samples (default 2).
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 1800
    - name: count
      description: >-
        Count reads in a consensus peakset for every sample (summarizeOverlaps,
        Union mode) and normalize for library size (CPM by default, RPKM also
        divides by peak width). Builds the consensus from peaks + indexed BAMs
        first, and populates a DiffBind DBA object when a sample sheet is present
        so differential binding can run next. Reports the count matrix shape and
        per-sample library sizes. The differential engines do their own internal
        normalization; this matrix is for exploration.
      run: [Rscript, ./scripts/count_reads.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A peaks directory/file plus BAMs, or a DiffBind sample sheet CSV with Peaks and bamReads.
          output:
            type: string
            format: project-path
            description: Where to save the counted DBA object / consensus and count tables as .rds.
          figures_dir:
            type: string
            format: project-path
            description: Directory for figures.
          peaks:
            type: array
            items:
              type: string
              format: input-path
            description: Explicit peak files; empty derives them from the sample sheet or input.
          bam:
            type: array
            items:
              type: string
              format: input-path
            description: Explicit indexed BAM files; empty derives them from the sample sheet, bam_dir or input.
          bam_dir:
            type: string
            format: input-path
            description: Directory of .bam files (with .bai indexes).
          sample_sheet:
            type: string
            format: input-path
            description: DiffBind-format CSV; null uses input when it is a sample sheet.
          paired_end:
            type: boolean
            description: Sequencing is paired-end (default false); affects counting.
          blacklist:
            type: string
            format: input-path
            description: User ENCODE blacklist BED; overrides the inlined hg38/mm10 lists.
          genome:
            type: string
            description: Genome build (hg38 or mm10) selecting the inlined ENCODE v2 blacklist.
          min_occupancy:
            type: integer
            minimum: 1
            description: Consensus occupancy threshold (default 2).
          normalize_method:
            type: string
            enum: [libsize, rpkm, none]
            description: Library-size CPM (default), RPKM (also divides by peak width), or none.
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 1800
    - name: diff
      description: >-
        Differential binding across conditions. method deseq2 (default) or edger
        runs DiffBind on a consensus count matrix built from a sample sheet plus
        peaks and BAMs; method csaw runs window-based analysis for broad marks
        (H3K27me3/H3K36me3/H3K9me3) whose signal spreads over kilobases. A contrast
        compares two conditions; with no explicit contrast every condition pair
        with enough replicates is tested. Reports the number of differential
        peaks/regions at the FDR and log2FC cutoffs, the top table and an MA plot.
        csaw is a gap on darwin-arm64 (no osx-arm64 bioconda build) — use DiffBind
        there. Start from a DiffBind sample sheet.
      run: [Rscript, ./scripts/diff_binding.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A DiffBind sample sheet CSV (SampleID/Condition/Replicate/Peaks/PeakCaller/bamReads).
          output:
            type: string
            format: project-path
            description: Where to save the analyzed DBA object and differential tables as .rds.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the MA plot.
          method:
            type: string
            enum: [deseq2, edger, csaw]
            description: Differential engine (default deseq2). csaw is window-based for broad marks; absent on darwin-arm64.
          sample_sheet:
            type: string
            format: input-path
            description: DiffBind-format CSV; null uses input.
          peaks:
            type: array
            items:
              type: string
              format: input-path
            description: Explicit peak files (DiffBind route); empty derives them from the sample sheet.
          bam:
            type: array
            items:
              type: string
              format: input-path
            description: Explicit indexed BAM files; empty derives them from the sample sheet, bam_dir or input.
          bam_dir:
            type: string
            format: input-path
            description: Directory of .bam files (with .bai indexes).
          paired_end:
            type: boolean
            description: Sequencing is paired-end (default false); affects counting and csaw fragment handling.
          blacklist:
            type: string
            format: input-path
            description: User ENCODE blacklist BED (DiffBind route); overrides the inlined hg38/mm10 lists.
          genome:
            type: string
            description: Genome build (hg38 or mm10) selecting the inlined ENCODE v2 blacklist.
          min_occupancy:
            type: integer
            minimum: 1
            description: Consensus occupancy threshold for the DiffBind route (default 2).
          contrast:
            type: array
            items:
              type: string
            description: Explicit contrast c(column, treated, control); empty tests all condition pairs with enough replicates.
          min_members:
            type: integer
            minimum: 1
            description: Minimum replicates per group for a contrast (default 2).
          fdr:
            type: number
            description: FDR cutoff for calling a peak differential (default 0.05).
          lfc:
            type: number
            description: Minimum absolute log2 fold change (default 0).
          top_n:
            type: integer
            minimum: 1
            description: Top differential peaks/regions printed and written (default 20).
          csaw_width:
            type: integer
            minimum: 1
            description: csaw window width in bp (default 150); larger for broad marks.
          csaw_spacing:
            type: integer
            minimum: 1
            description: csaw window spacing in bp (default 50); also the region-merge tolerance.
          csaw_filter:
            type: string
            description: csaw filterWindows type (default global); global background filtering suits broad marks.
          fragment_size:
            type: integer
            minimum: 1
            description: csaw single-end fragment extension in bp (default 200); ignored for paired-end.
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: annotate
      description: >-
        Annotate peaks to genomic features with ChIPseeker (promoter/exon/intron/
        intergenic distribution, nearest gene, distance to TSS) and draw a TSS
        coverage profile and peak heatmap. The lean annotation strategy means the
        target-genome TxDb comes from the user: a GTF/GFF3 (gtf) or a prebuilt TxDb
        .sqlite (txdb); without one this step reports a gap instead of guessing.
        Gene symbols need a user org.db; without it annotation reports gene ids.
        Accepts a peak file/dir, a sample sheet, or a saved .rds (DBA or GRanges).
      run: [Rscript, ./scripts/annotate_peaks.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A peaks file/dir, a sample sheet, or a saved .rds (DBA or GRanges) from an earlier step.
          output:
            type: string
            format: project-path
            description: Where to save the annotated object and annotation table as .rds.
          figures_dir:
            type: string
            format: project-path
            description: Directory for feature-distribution, TSS-distance, TSS-profile and peak-heatmap figures.
          peaks:
            type: array
            items:
              type: string
              format: input-path
            description: Explicit peak files; empty derives them from the sample sheet or input.
          sample_sheet:
            type: string
            format: input-path
            description: DiffBind-format CSV; null uses input when it is a sample sheet.
          gtf:
            type: string
            format: input-path
            description: User GTF/GFF3 to build a TxDb (makeTxDbFromGFF); required if txdb is absent.
          txdb:
            type: string
            format: input-path
            description: User prebuilt TxDb .sqlite (loadDb); overrides gtf.
          org_db:
            type: string
            description: Installed org.db package name (e.g. org.Hs.eg.db) for gene symbols; absent reports gene ids.
          tss_upstream:
            type: integer
            minimum: 0
            description: TSS region upstream extent in bp (default 3000).
          tss_downstream:
            type: integer
            minimum: 0
            description: TSS region downstream extent in bp (default 3000).
          blacklist:
            type: string
            format: input-path
            description: User ENCODE blacklist BED; overrides the inlined hg38/mm10 lists.
          genome:
            type: string
            description: Genome build (hg38 or mm10) selecting the inlined ENCODE v2 blacklist.
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: visualize
      description: >-
        Exploratory views of a peak x sample count matrix: sample PCA on
        log-normalized counts (the first check of whether the dominant variation is
        condition or a nuisance) and a peak-overlap UpSet plot. Entry point is
        either a count table (counts or a .csv/.tsv input) or peaks + BAMs that get
        counted here. Reports variance explained per PC and the sample scores.
        Start from a count table, or peaks + BAMs.
      run: [Rscript, ./scripts/visualize.R]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A peak x sample count table (.csv/.tsv), or peaks + BAMs / a sample sheet to count here.
          output:
            type: string
            format: project-path
            description: Where to save the object and normalized count table as .rds.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the sample-PCA and peak-UpSet figures.
          counts:
            type: string
            format: input-path
            description: Peak x sample count table (first column peak id, header row of samples); null uses input when it is a table.
          peaks:
            type: array
            items:
              type: string
              format: input-path
            description: Explicit peak files for the counting route; empty derives them from the sample sheet or input.
          bam:
            type: array
            items:
              type: string
              format: input-path
            description: Explicit indexed BAM files for the counting route.
          bam_dir:
            type: string
            format: input-path
            description: Directory of .bam files (with .bai indexes).
          sample_sheet:
            type: string
            format: input-path
            description: DiffBind-format CSV; null uses input when it is a sample sheet.
          paired_end:
            type: boolean
            description: Sequencing is paired-end (default false); affects counting.
          blacklist:
            type: string
            format: input-path
            description: User ENCODE blacklist BED; overrides the inlined hg38/mm10 lists.
          genome:
            type: string
            description: Genome build (hg38 or mm10) selecting the inlined ENCODE v2 blacklist.
          min_occupancy:
            type: integer
            minimum: 1
            description: Consensus occupancy threshold for the counting route (default 2).
          normalize_method:
            type: string
            enum: [libsize, rpkm, none]
            description: Library-size CPM (default), RPKM, or none.
          sample_column:
            type: string
            description: Sample sheet column naming each sample, for PCA colouring (default SampleID).
          seed:
            type: integer
            description: Random seed (default 42).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 1800
    - name: notebook
      description: >-
        Render a runnable R analysis notebook into the project from the fragment
        library, then open it with Phi's notebook tools and run cells on demand.
        The template is rendered, not copied: placeholders are filled from the
        input data and the requested steps, and every tunable value lands in one
        config cell at the top that the analysis cells reference. Pick a preset
        with module (1 TF/narrow DiffBind+DESeq2, 2 broad csaw, 3 CUT&Tag edgeR,
        4 consensus, 5 annotation, 6 count exploration) or name steps directly;
        dependencies are added automatically. inject adds an external R or Python
        script as a code cell; merge splices in another notebook, whole or a cell
        range. All sources must match the notebook language (one kernel per
        notebook). Starts from called peaks / counts / BAMs: raw-read alignment and
        peak calling are the wrapper tree's job, not a notebook cell. Use self_check
        to execute every cell in the output directory now and fail loudly if one
        breaks.
      run: [python, ./scripts/build_notebook.py]
      args:
        type: object
        additionalProperties: false
        required: [output, data]
        properties:
          module:
            type: integer
            enum: [1, 2, 3, 4, 5, 6]
            description: Analysis preset. Pass module or steps, not both.
          steps:
            type: array
            items:
              type: string
              enum: [setup, load-peaks, load-bam, load-counts, load-sample-sheet, qc-frip, qc-reproducibility, blacklist-filter, consensus, count, normalize, pca, diff-deseq2, diff-edger, csaw, annotate, coverage-profile, upset, save]
            description: >-
              Steps to render instead of a preset; declared dependencies are added
              automatically. Pass module or steps, not both.
          output:
            type: string
            format: project-path
            description: Path of the .ipynb to write inside the project.
          data:
            type: string
            format: input-path
            description: Called peaks (file/dir), a DiffBind sample sheet, a count table, BAMs, or a saved .rds the notebook starts from.
          language:
            type: string
            enum: [r, python]
            description: >-
              Notebook language and kernel (default r). The R fragment library backs
              module/steps; a python notebook is assembled from inject and merge
              sources only.
          inject:
            type: array
            items:
              type: string
              format: input-path
            description: R (.R) or Python (.py) scripts to add verbatim as code cells, in order.
          merge:
            type: array
            items:
              type: string
            description: >-
              Notebooks to splice in after the rendered steps, in order. Whole
              notebook (path.ipynb) or a cell range (path.ipynb[0:6], indices over
              all cells). Must match the notebook language.
          config:
            type: string
            format: input-path
            description: JSON file of parameter overrides; keys mirror the parameter names below.
          project:
            type: string
            description: Project name recorded on the object; defaults to the data folder name.
          output_object:
            type: string
            description: Where the notebook saves its final object, relative to the notebook folder (default analyzed.rds).
          sample_column:
            type: string
            description: Sample sheet column naming each sample (default SampleID).
          kernel:
            type: string
            description: Kernel name written to the notebook (default plugin-chromatin-profiling-chip / phi-python).
          self_check:
            type: boolean
            description: Execute every cell once in the output directory and fail if one breaks.
          seed:
            type: integer
            description: Random seed (default 42).
          figures_dir:
            type: string
            description: Directory for PNG figures (default figures); null disables saving.
          peak_paths:
            type: array
            items:
              type: string
            description: Explicit peak files; empty uses data_path (file or directory).
          count_path:
            type: string
            description: Peak x sample count table; null uses data_path when it is a table.
          bam_paths:
            type: array
            items:
              type: string
            description: Explicit indexed BAM files; empty uses bam_dir or data_path.
          bam_dir:
            type: string
            description: Directory of .bam files (with .bai indexes).
          sample_sheet_path:
            type: string
            description: DiffBind-format CSV; null uses data_path when it is a sample sheet.
          paired_end:
            type: boolean
            description: Sequencing is paired-end (default false).
          blacklist_path:
            type: string
            description: User ENCODE blacklist BED; overrides the inlined hg38/mm10 lists.
          genome:
            type: string
            description: Genome build (hg38 or mm10) selecting the inlined ENCODE v2 blacklist.
          min_occupancy:
            type: integer
            description: Consensus occupancy threshold (default 2).
          normalize_method:
            type: string
            description: libsize (CPM, default), rpkm, or none.
          contrast:
            type: array
            items:
              type: string
            description: Explicit contrast c(column, treated, control); empty tests all condition pairs.
          min_members:
            type: integer
            description: Minimum replicates per group for a DiffBind contrast (default 2).
          fdr_threshold:
            type: number
            description: FDR cutoff for a differential peak (default 0.05).
          lfc_threshold:
            type: number
            description: Minimum absolute log2 fold change (default 0).
          top_n:
            type: integer
            description: Top differential peaks/regions printed and written (default 20).
          csaw_width:
            type: integer
            description: csaw window width in bp (default 150).
          csaw_spacing:
            type: integer
            description: csaw window spacing in bp (default 50); also the region-merge tolerance.
          csaw_filter:
            type: string
            description: csaw filterWindows type (default global).
          fragment_size:
            type: integer
            description: csaw single-end fragment extension in bp (default 200).
          gtf_path:
            type: string
            description: User GTF/GFF3 to build a TxDb for annotation.
          txdb_path:
            type: string
            description: User prebuilt TxDb .sqlite; overrides gtf_path.
          org_db:
            type: string
            description: Installed org.db package name for gene symbols.
          tss_upstream:
            type: integer
            description: TSS region upstream extent in bp (default 3000).
          tss_downstream:
            type: integer
            description: TSS region downstream extent in bp (default 3000).
      approval: write
      output: ./schemas/notebook-result.json
      timeoutSeconds: 1800
    - name: env_probe
      description: >-
        Report which ChIP/CUT&Tag downstream R packages the current environment
        provides, list the packages each analysis preset needs, and say which
        preset cores can run here and which optional sections cannot. Also lists
        structural gaps (csaw on darwin-arm64, IDR, profileplyr, motif enrichment,
        bundled org.db, upstream peak calling) with a concrete next step. Call this
        before promising a workflow such as broad-mark csaw on Apple Silicon. The
        plugin environment is read-only and prebuilt, so a missing package needs a
        project environment (env_request), not an in-environment install.
      run: [Rscript, ./scripts/env_probe.R]
      args:
        type: object
        additionalProperties: false
        properties:
          module:
            type: array
            items:
              type: integer
            description: Restrict the report to these preset numbers (1-6).
      approval: read
      output: ./schemas/env-probe.json
      timeoutSeconds: 300
---

# Chromatin profiling (bulk ChIP / CUT&Tag) workflow

Run a reproducible downstream analysis on called peaks, one experiment at a time.
Every step reads peaks / counts / BAMs (or a saved `.rds`) and writes reviewable
tables and figures, so steps chain, can be re-run in isolation, and leave an audit
trail. The R fragment library (`assets/fragments/*.R`) is the single analysis
source: the script tools and the rendered notebook run the same fragments.

## Before the first step

- Call `chip_inspect` on an unfamiliar input. It reports the input kind, shape,
  per-sample peak counts, sample-sheet columns/conditions, BAM index status, and
  which steps the input is ready for.
- Call `chip_env_probe` before promising csaw on Apple Silicon, IDR, motif
  enrichment or gene symbols. It states what the locked environment provides and
  what is a gap, with a next step. Never install into the plugin environment.
- This plugin is **pure downstream**. It starts from called peaks
  (BED/narrowPeak/broadPeak), a DiffBind sample sheet, a count table, aligned BAMs
  or a bigWig. FASTQ → alignment → peak calling belongs to the nf-core/chipseq
  wrapper (ChIP) or the nf-core/seacr/callpeak module wrapper (CUT&Tag).

## The chain

1. **QC** (`chip_qc`) — FRiP per sample and replicate peak-overlap reproducibility.
   Read the FRiP distribution before judging enrichment; TF ChIP ideally > 1%,
   broad marks lower, CUT&Tag often much higher. IDR is the formal reproducibility
   test and is a gap; overlap fraction is a transparent proxy.
2. **Consensus** (`chip_consensus`) — blacklist filter, then a merged peak universe
   present in ≥ `min_occupancy` samples. This is the region set counting tests over.
3. **Count** (`chip_count`) — reads per consensus peak per sample, library-size
   normalized for exploration. Populates a DiffBind DBA when a sample sheet exists.
4. **Differential binding** (`chip_diff`) — `method deseq2|edger` (DiffBind, TF /
   narrow marks and CUT&Tag) or `method csaw` (window-based, broad marks). Read the
   MA plot and the FDR/log2FC counts; choose the contrast deliberately.
5. **Annotate** (`chip_annotate`) — ChIPseeker feature distribution, TSS distance,
   coverage profile and heatmap against a **user** GTF/GFF3 or TxDb. Lean annotation:
   no target-genome TxDb is bundled; without one this step reports a gap.
6. **Visualize** (`chip_visualize`) — sample PCA and overlap UpSet on a count matrix.

`chip_consensus`, `chip_count` and `chip_diff` accept a DiffBind sample sheet as
`input`; the sheet's `Peaks`/`bamReads` columns drive peak and BAM loading. Relative
paths in the sheet resolve against the sheet's directory.

## Rendering a notebook

`chip_notebook` renders a runnable R notebook into the project (never copies a
template directory). Pick a preset with `module` or name `steps`; dependencies are
added automatically. Every tunable value lands in one config cell. The inlined
ENCODE blacklist travels with the `blacklist-filter` fragment, so the notebook
filters offline in the user's directory. Use `self_check` to execute every cell once
and fail loudly if one breaks. The kernel is `plugin-chromatin-profiling-chip`.

## Rules

- **Pure downstream.** Never expect FASTQ or do alignment/peak calling here. Point
  the user to the nf-core/chipseq wrapper (ChIP) or their SEACR run (CUT&Tag).
- **Read-only environment.** A missing package or genome resource is a fact to
  report with a next step (project environment via `env_request`, or a user
  `input-path`), never an in-environment install or large download.
- **Lean annotation.** The target-genome TxDb comes from the user (`gtf`/`txdb`).
  The hg19 TxDb that ships as a ChIPseeker conda build dependency is **never** the
  default annotation source for the target genome. Gene symbols need a user org.db.
- **Blacklist is genome-build specific.** hg38 regions are wrong on mm10. Pass
  `genome` (hg38/mm10) for the inlined ENCODE v2 lists, or `blacklist` for another
  build. With neither, filtering skips and says so.
- **DiffBind vs csaw is a biology decision.** TF/narrow marks and CUT&Tag → DiffBind
  (DESeq2/edgeR). Broad marks (H3K27me3/H3K36me3/H3K9me3) → csaw windows. csaw is
  absent on darwin-arm64; use DiffBind there or run on linux-x64/darwin-x64.
- **Replicates are the experimental unit.** Fewer than `min_members` (default 2)
  biological replicates per condition → no condition-level claim. Normalized/exploratory
  counts never feed the differential test; the engines normalize internally.
- **The JSON return is the audit record.** Quote peak counts, FRiP, occupancy,
  library sizes, differential counts and what could not run; do not recompute them.

## Reading the results

Each tool emits one JSON object on stdout (diagnostics on stderr). `chip_qc` gives
`frip` and `reproducibility`; `chip_consensus` gives `n_consensus` and
`occupancy_per_sample`; `chip_count` gives `library_sizes`; `chip_diff` gives
`method`, `n_significant`, `n_tested` and the top `table` (or an `error` for the csaw
platform gap); `chip_annotate` gives `feature_distribution` and `annotation_source`;
`chip_visualize` gives `pca.variance_explained` and `scores`. Figures land in
`figures_dir`; tables (consensus.bed, counts.csv, diff_binding.csv, csaw_regions.csv,
annotation.csv) land beside the saved object. See `references/analysis-playbook.md`
for the call → read → decide → red-line loop and `references/step-notes/` for each
step's purpose and pitfalls.