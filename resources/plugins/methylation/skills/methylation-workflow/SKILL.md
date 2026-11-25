---
name: methylation-workflow
description: >-
  Downstream DNA methylation workflow in Python for WGBS, RRBS, EM-seq
  (nf-core/methylseq, Bismark .bismark.cov.gz, MethylDackel .bedGraph), and
  Illumina 450K/EPIC beta-value matrices: inspect inputs and coverage; merge
  per-sample cytosine coverage files into depth-filtered beta and coverage
  matrices; run global beta bimodality and logit M-value PCA QC; call
  confounder-adjusted empirical-Bayes Differentially Methylated Cytosines (DMCs)
  and sliding-window Differentially Methylated Regions (DMRs, Stouffer regional
  FDR + delta-beta cutoff); and annotate against GTF promoters/gene bodies and
  CpG islands/shores/shelves with GMT ORA. Use for DNA甲基化, methylation, WGBS,
  RRBS, EM-seq, Bismark, 450K, EPIC, beta值, M-value, DMC, DMR, 差异甲基化, or CpG岛.
phi:
  environment: plugin:methylation
  attachTo:
    - Methylation
  scripts:
    - name: inspect
      description: >-
        Inspect per-sample Bismark .cov(.gz) / MethylDackel .bedGraph(.gz) files
        or a CpG/probe x sample methylation matrix (beta or M-value) and optional
        sample metadata: report input format, CpG count, value range, coverage
        summary, and metadata sample alignment. Call first on unfamiliar inputs.
      run: [python, ./scripts/meth_inspect.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A methylation matrix TSV/CSV, a Bismark/bedGraph file, or a directory of per-sample coverage files.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV.
          sample_col:
            type: string
            description: Sample ID column in metadata (default first column).
          group_col:
            type: string
            description: Optional condition column in metadata.
      approval: read
      output: ./schemas/inspect-result.json
      timeoutSeconds: 180
    - name: merge
      description: >-
        Merge per-sample Bismark (.bismark.cov / .cov.gz) or MethylDackel
        (.bedGraph / .bedGraph.gz) cytosine methylation reports into a unified
        chr:pos x sample beta-value matrix and read-coverage matrix, filtering by
        minimum per-sample read depth and minimum sample call rate.
      run: [python, ./scripts/meth_merge.py]
      args:
        type: object
        additionalProperties: false
        required: [inputs, output_dir]
        properties:
          inputs:
            type: array
            items:
              type: string
            description: Per-sample coverage files (optionally sample=path) or a single directory path.
          min_cov:
            type: integer
            minimum: 1
            description: Minimum total reads (meth + unmeth) required to call a CpG in a sample (default 5).
          min_samples_frac:
            type: number
            description: Minimum fraction of samples with valid coverage >= min_cov to retain a CpG (default 0.7).
          output_dir:
            type: string
            format: project-path
            description: Output directory for beta_matrix.tsv and coverage_matrix.tsv.
      approval: write
      output: ./schemas/merge-result.json
      timeoutSeconds: 600
    - name: qc
      description: >-
        Run global DNA methylation QC on a beta or M-value matrix: summarize
        per-sample mean/median beta, hypo- (<0.2), intermediate (0.2-0.8), and
        hyper-methylated (>0.8) fractions, missingness, logit M-value PCA, and
        outlier flags, writing sample_qc.tsv and meth-qc.png.
      run: [python, ./scripts/meth_qc.py]
      args:
        type: object
        additionalProperties: false
        required: [matrix, output_dir]
        properties:
          matrix:
            type: string
            format: input-path
            description: CpG/probe x sample beta-value (0..1) or M-value matrix TSV/CSV.
          coverage_matrix:
            type: string
            format: input-path
            description: Optional CpG x sample read-depth matrix TSV/CSV from meth_merge.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV for PCA group coloring.
          sample_col:
            type: string
            description: Sample ID column in metadata (default first column).
          group_col:
            type: string
            description: Optional group column in metadata.
          output_dir:
            type: string
            format: project-path
            description: Output directory for sample_qc.tsv and meth-qc.png.
      approval: write
      output: ./schemas/qc-result.json
      timeoutSeconds: 300
    - name: dmr
      description: >-
        Call Differentially Methylated Cytosines (DMCs) via empirical-Bayes
        moderated linear regression on logit M-values (supporting confounder
        covariates) paired with exact biological delta-beta thresholds, and
        merge spatially adjacent directional CpGs into Differentially Methylated
        Regions (DMRs) with Stouffer's Z combined p-values and BH FDR.
      run: [python, ./scripts/meth_dmr.py]
      args:
        type: object
        additionalProperties: false
        required: [matrix, metadata, group_col, case, control, output_dir]
        properties:
          matrix:
            type: string
            format: input-path
            description: CpG/probe x sample beta-value (0..1) or M-value matrix TSV/CSV.
          metadata:
            type: string
            format: input-path
            description: Sample metadata TSV/CSV.
          sample_col:
            type: string
            description: Sample ID column in metadata (default first column).
          group_col:
            type: string
            description: Grouping column in metadata.
          case:
            type: string
            description: Case group value (delta_beta = mean_case - mean_control).
          control:
            type: string
            description: Control/reference group value.
          covariates:
            type: array
            items:
              type: string
            description: Optional confounder columns in metadata (e.g. age, sex, batch, purity).
          probe_coords:
            type: string
            format: input-path
            description: Optional probe coordinate TSV/BED (probe_id, chrom, pos) when matrix IDs are cg* probes rather than chr:pos.
          delta_beta:
            type: number
            description: Minimum absolute delta-beta threshold for biological significance (default 0.10).
          fdr_cutoff:
            type: number
            description: BH FDR cutoff for DMCs and DMRs (default 0.05).
          max_dist:
            type: integer
            minimum: 50
            description: Maximum base-pair gap to merge adjacent directional CpGs into one DMR (default 500).
          min_cpgs:
            type: integer
            minimum: 2
            description: Minimum CpGs required per called DMR (default 3).
          output_dir:
            type: string
            format: project-path
            description: Output directory for dmc_results.tsv, dmr_regions.tsv, and meth-dmr.png.
      approval: write
      output: ./schemas/dmr-result.json
      timeoutSeconds: 600
    - name: annotate
      description: >-
        Annotate DMRs or DMCs against a reference GTF/GFF3 (Promoter, GeneBody,
        Downstream, Distal_Intergenic + nearest gene symbol & signed TSS
        distance) and optional CpG Island BED (Island, Shore +-2kb, Shelf 2-4kb,
        OpenSea), with optional GMT hypergeometric ORA against the assayed gene
        universe.
      run: [python, ./scripts/meth_annotate.py]
      args:
        type: object
        additionalProperties: false
        required: [input, gtf, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: DMR table (dmr_regions.tsv) or coordinate-bearing DMC table (dmc_results.tsv).
          gtf:
            type: string
            format: input-path
            description: Reference GTF/GFF3 (.gtf/.gtf.gz/.gff3) supplying gene TSS and body coordinates.
          cpg_islands:
            type: string
            format: input-path
            description: Optional CpG Island BED file (chrom, start, end) for Island/Shore/Shelf classification.
          gmt:
            type: string
            format: input-path
            description: Optional GMT gene-set file for ORA on hyper/hypomethylated target genes.
          all_tested_table:
            type: string
            format: input-path
            description: Optional full dmc_results.tsv to define the exact assayed CpG gene universe for ORA.
          promoter_upstream:
            type: integer
            description: Upstream promoter window in bp from TSS (default 2000).
          promoter_downstream:
            type: integer
            description: Downstream promoter window in bp from TSS (default 500).
          output_dir:
            type: string
            format: project-path
            description: Output directory for annotated table, optional ORA table, and genomic context plot.
      approval: write
      output: ./schemas/annotate-result.json
      timeoutSeconds: 600
---

# DNA Methylation Downstream Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) before
choosing coverage cutoffs, $\Delta\beta$ thresholds, confounder adjustment, or
interpreting promoter vs gene-body methylation.

## Chain

```text
FASTQ ──(nf-core/methylseq / bismark / methyldackel)──▶ per-sample .cov.gz / .bedGraph
  ──meth_inspect──▶ verify format, depth, beta range [0,1] & metadata alignment
  ──meth_merge────▶ coverage-filtered beta_matrix.tsv + coverage_matrix.tsv
  ──meth_qc───────▶ bimodal beta distribution & logit M-value PCA
  ──meth_dmr──────▶ moderated M-value linear model + |Δβ| DMC & Stouffer DMRs
  ──meth_annotate─▶ Promoter/GeneBody + CpG Island/Shore/Shelf + assayed-universe ORA
```

## Hard Red Lines

1. **Statistics on $M$-values, effect sizes on $\Delta\beta$**: never run linear regression directly on bounded $\beta \in [0, 1]$ near $0/1$, and never interpret tiny $|\Delta\beta| < 0.02$ as major biological reprogramming without cell-purity control.
2. **Enforce read-depth floors on WGBS/RRBS**: a CpG with $2\times$ coverage can only take $\beta \in \{0, 0.5, 1.0\}$; require `--min_cov >= 5` (`>= 10` preferred when depth allows).
3. **Promoter hypermethylation $\neq$ gene-body hypermethylation**: CpG-island promoter hypermethylation typically represses transcription, whereas gene-body methylation often correlates positively with active transcription.
