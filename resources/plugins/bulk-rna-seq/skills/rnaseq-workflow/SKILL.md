---
name: rnaseq-workflow
description: >-
  Downstream bulk RNA-seq interpretation and pathway analysis: inspect count
  or normalized matrices, DE tables, sample metadata, and GMT gene sets
  (detecting TPM-as-count traps and batch-condition collinearity); run sample
  QC/PCA/correlation outlier diagnostics; integrate multi-contrast DESeq2 /
  edgeR / limma DEG tables; perform offline hypergeometric ORA with explicit
  background universe; run deterministic pre-ranked GSEA (ES/NES/FDR +
  leading-edge + enrichment curves); and compute single-sample pathway
  activity scores. Use for bulk RNA-seq, 转录组, 差异表达下游, GSEA, ORA,
  通路富集, 多对比差异整合, or pathway activity scoring. Raw FASTQ alignment
  runs in nf-core/rna-seq and core DE models run in
  local/differential-expression wrappers.
phi:
  environment: plugin:rnaseq
  attachTo:
    - BulkRnaSeq
  scripts:
    - name: inspect
      description: >-
        Inspect an expression matrix (inferring raw_counts vs tpm_or_cpm vs
        log_normalized), sample metadata (sample alignment, replicate counts,
        batch-group collinearity), existing DE table, and/or GMT gene-set file.
        Call first before launching DESeq2/edgeR or pathway analysis.
      run: [python, ./scripts/rna_inspect.py]
      args:
        type: object
        additionalProperties: false
        properties:
          matrix:
            type: string
            format: input-path
            description: Gene x sample count or normalized expression matrix TSV/CSV.
          metadata:
            type: string
            format: input-path
            description: Sample metadata TSV/CSV.
          sample_col:
            type: string
            description: Optional metadata column holding sample IDs (default first column).
          group_col:
            type: string
            description: Primary condition column in metadata (default condition).
          batch_col:
            type: string
            description: Optional batch column to test for collinearity with group_col.
          de_table:
            type: string
            format: input-path
            description: Optional DESeq2/edgeR/limma result table.
          gmt:
            type: string
            format: input-path
            description: Optional GMT gene-set file.
      approval: read
      output: ./schemas/inspect-result.json
      timeoutSeconds: 300
    - name: qc
      description: >-
        Run sample-level bulk RNA-seq QC on a count or normalized expression
        matrix: library sizes, detected gene counts, log2(CPM+1) PCA on top
        variable genes, sample-to-sample Pearson correlation heatmap, and
        outlier flagging.
      run: [python, ./scripts/rna_qc.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Gene x sample expression matrix TSV/CSV.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV.
          sample_col:
            type: string
            description: Optional metadata sample ID column.
          group_col:
            type: string
            description: Metadata column used to color PCA points (default condition).
          top_var:
            type: integer
            minimum: 10
            description: Number of top variable genes for PCA (default 500).
          min_corr_outlier:
            type: number
            description: Mean Pearson correlation floor below which a sample is flagged (default 0.80).
          output_dir:
            type: string
            format: project-path
            description: Directory for QC tables and PCA/correlation figures.
      approval: write
      output: ./schemas/qc-result.json
      timeoutSeconds: 600
    - name: compare
      description: >-
        Integrate two or more DESeq2/edgeR/limma contrast tables into one
        multi-contrast DEG matrix with direction-concordance labels
        (shared_up, shared_down, discordant, contrast_specific), pairwise
        Jaccard overlap and log2FC correlation table, and comparison figure.
      run: [python, ./scripts/rna_compare.py]
      args:
        type: object
        additionalProperties: false
        required: [inputs, output_dir]
        properties:
          inputs:
            type: array
            items:
              type: string
            description: Two or more DE tables, optionally formatted as label=path.
          padj:
            type: number
            description: Adjusted p-value threshold (default 0.05).
          log2fc:
            type: number
            description: Minimum absolute log2 fold-change threshold (default 1.0).
          output_dir:
            type: string
            format: project-path
            description: Directory for integrated DEG tables and comparison plot.
      approval: write
      output: ./schemas/compare-result.json
      timeoutSeconds: 600
    - name: ora
      description: >-
        Offline hypergeometric Over-Representation Analysis (ORA) from a full
        DE table (automatically supplying the expressed background universe) or
        explicit --genes + --background lists against a GMT file, with BH FDR
        and enrichment dotplot.
      run: [python, ./scripts/rna_ora.py]
      args:
        type: object
        additionalProperties: false
        required: [gmt, output_dir]
        properties:
          de_table:
            type: string
            format: input-path
            description: Full DESeq2/edgeR/limma table (supplies foreground DEGs and universe).
          genes:
            type: string
            format: input-path
            description: Foreground gene list (one per line) when de_table is omitted.
          background:
            type: string
            format: input-path
            description: Explicit background universe gene list (mandatory with --genes).
          gmt:
            type: string
            format: input-path
            description: GMT gene-set file (e.g., MSigDB Hallmark/GO/KEGG/Reactome).
          direction:
            type: string
            enum: [both, up, down]
            description: DEG direction to test from de_table (default both).
          padj:
            type: number
            description: DE padj cutoff (default 0.05).
          log2fc:
            type: number
            description: DE |log2FC| cutoff (default 1.0).
          min_set_size:
            type: integer
            minimum: 1
            description: Minimum gene-set size overlapping universe (default 5).
          max_set_size:
            type: integer
            minimum: 5
            description: Maximum gene-set size overlapping universe (default 500).
          fdr:
            type: number
            description: Pathway BH FDR cutoff (default 0.05).
          output_dir:
            type: string
            format: project-path
            description: Directory for ORA enrichment table and dotplot.
      approval: write
      output: ./schemas/ora-result.json
      timeoutSeconds: 600
    - name: gsea
      description: >-
        Offline deterministic pre-ranked GSEA (weighted Kolmogorov-Smirnov ES,
        seeded gene-set permutation NES, p-value, BH FDR, leading-edge genes,
        and running-sum enrichment curves) from a full unfiltered DE table and
        a GMT file.
      run: [python, ./scripts/rna_gsea.py]
      args:
        type: object
        additionalProperties: false
        required: [de_table, gmt, output_dir]
        properties:
          de_table:
            type: string
            format: input-path
            description: Full unfiltered DE table (DESeq2/edgeR/limma).
          gmt:
            type: string
            format: input-path
            description: GMT gene-set file.
          rank_by:
            type: string
            enum: [stat, log2fc]
            description: Ranking metric column (default stat; falls back to signed -log10(p) if absent).
          weight:
            type: number
            description: GSEA weighting exponent p (default 1.0).
          min_set_size:
            type: integer
            minimum: 1
            description: Minimum gene-set size overlapping ranked list (default 5).
          max_set_size:
            type: integer
            minimum: 5
            description: Maximum gene-set size overlapping ranked list (default 500).
          permutations:
            type: integer
            minimum: 100
            description: Number of gene-set permutations (default 500).
          fdr:
            type: number
            description: Pathway BH FDR cutoff (default 0.05).
          seed:
            type: integer
            description: RNG seed for permutations (default 12345).
          output_dir:
            type: string
            format: project-path
            description: Directory for GSEA results table and enrichment curves.
      approval: write
      output: ./schemas/gsea-result.json
      timeoutSeconds: 1200
    - name: pathway_score
      description: >-
        Compute single-sample pathway activity scores (mean per-gene z-score
        across samples for each GMT pathway), optional group comparison
        (Mann-Whitney U / Kruskal-Wallis + BH FDR), and a pathway-by-sample
        activity heatmap.
      run: [python, ./scripts/rna_pathway_score.py]
      args:
        type: object
        additionalProperties: false
        required: [input, gmt, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Gene x sample expression matrix TSV/CSV (raw counts or normalized).
          gmt:
            type: string
            format: input-path
            description: GMT gene-set file.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata for group testing.
          sample_col:
            type: string
            description: Optional metadata sample ID column.
          group_col:
            type: string
            description: Metadata column defining groups (default condition).
          min_set_size:
            type: integer
            minimum: 1
            description: Minimum matched genes required to score a pathway (default 3).
          output_dir:
            type: string
            format: project-path
            description: Directory for score matrix, group tests, and heatmap.
      approval: write
      output: ./schemas/pathway-score-result.json
      timeoutSeconds: 600
---

# Bulk RNA-seq — downstream interpretation & pathway chain

Read [references/analysis-playbook.md](references/analysis-playbook.md) before
choosing between `DESeq2`, `edgeR`, `limma`, `timeseries`, `rna_ora`, and
`rna_gsea`.

## Chain

```
FASTQ ──(nf-core/rna-seq wrapper)──▶ gene counts
counts + metadata ──rna_inspect──▶ verify raw counts & batch-group design
                  ──rna_qc───────▶ library sizes, PCA, correlation & outlier check
                  ──(local/differential-expression: deseq2 / edger / limma / timeseries)──▶ DE table(s)
multiple DE tables ──rna_compare──▶ shared / discordant / specific DEGs
full DE table + GMT ──┬─ rna_ora  (hypergeometric ORA with explicit universe)
                      └─ rna_gsea (threshold-free pre-ranked GSEA + leading-edge)
expression + GMT ──rna_pathway_score──▶ sample-level pathway activity matrix + heatmap
```

## Hard red lines

1. **Never pass TPM/FPKM/VST to DESeq2 or edgeR**: verify `inferred_kind == "raw_counts"` via `rna_inspect` first.
2. **Check batch-condition collinearity**: if `rna_inspect` flags `confounded: true`, do not silently fit `~ batch + condition` (rank-deficient).
3. **Explicit background universe in ORA**: pass the full DE table (or explicit `--background` of expressed genes), never an unstated genome-wide default.
4. **Unfiltered ranked table in GSEA**: never pre-filter the DE table to `padj < 0.05` before calling `rna_gsea`.
5. **Replicate floor**: warn explicitly whenever any group has `n < 3` replicates.
