---
name: rnaseq-workflow
description: >-
  Downstream bulk RNA-seq interpretation, identifier/GEO bridging, co-expression
  networks, TME deconvolution, and pathway analysis: import NCBI GEO Series
  Matrix cohorts; strip Ensembl versions and map gene IDs via GTF/GFF3 or tables;
  inspect matrices, DE tables, metadata, and GMT sets; run sample QC/PCA; compare
  multi-contrast DEGs; run hypergeometric ORA, pre-ranked GSEA, single-sample
  pathway scoring, ULM/MLM TF and pathway regulon inference, rMATS/SUPPA2
  differential alternative splicing harmonization, WGCNA-style weighted gene
  co-expression module & hub-gene detection, and ESTIMATE/NNLS tumor
  microenvironment immune deconvolution. Use for bulk RNA-seq, 转录组,
  基因ID转换, GEO数据导入, 差异表达下游, GSEA, ORA, regulon, 转录因子活性,
  可变剪接, rMATS, SUPPA2, WGCNA, 共表达网络, 枢纽基因, 免疫浸润,
  肿瘤微环境去卷积. Raw alignment runs on nf-core/rna-seq and core DE on
  local/differential-expression wrappers.
phi:
  environment: plugin:rnaseq
  attachTo:
    - BulkRnaSeq
  scripts:
    - name: id_map
      description: >-
        Clean versioned Ensembl/RefSeq gene IDs (e.g. ENSG00000141510.18 ->
        ENSG00000141510), optionally filter by GTF biotype (e.g. protein_coding),
        and map gene identifiers using a reference GTF/GFF3 or custom 2-column
        mapping table with collision-aware aggregation (sum for raw counts,
        max_mean for normalized matrices, min_pvalue for DE tables).
      run: [python, ./scripts/rna_id_map.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output]
        properties:
          input:
            type: string
            format: input-path
            description: Expression matrix or DE result table TSV/CSV.
          output:
            type: string
            format: project-path
            description: Output mapped TSV/CSV path.
          gtf:
            type: string
            format: input-path
            description: Optional reference GTF/GFF3 (.gtf/.gtf.gz/.gff3) for gene_id -> gene_name mapping.
          mapping_table:
            type: string
            format: input-path
            description: Optional 2-column TSV/CSV mapping file.
          from_col:
            type: string
            description: Optional source column name in mapping_table.
          to_col:
            type: string
            description: Optional target column name in mapping_table.
          id_col:
            type: string
            description: Optional gene ID column in input (default first column).
          keep_version:
            type: boolean
            description: Retain Ensembl/RefSeq '.version' suffixes instead of stripping them before lookup (default false).
          table_kind:
            type: string
            enum: [auto, matrix, de]
            description: Input table type (default auto).
          agg_method:
            type: string
            enum: [auto, sum, mean, max_mean, min_pvalue, first]
            description: Many-to-one symbol collision aggregation strategy (default auto).
          unmapped:
            type: string
            enum: [keep, drop]
            description: Whether to retain cleaned unmapped IDs or drop them (default keep).
          biotype:
            type: string
            description: Optional GTF gene_biotype filter (e.g. protein_coding).
      approval: write
      output: ./schemas/id-map-result.json
      timeoutSeconds: 300
    - name: geo_import
      description: >-
        Parse an NCBI GEO Series Matrix file (GSE*_series_matrix.txt or .txt.gz)
        into a structured samplesheet.tsv (parsing !Sample_characteristics_ch1
        key:value pairs into clean columns) and expression_matrix.tsv, with
        optional GPL platform probe-to-symbol mapping and value-scale detection.
      run: [python, ./scripts/rna_geo_import.py]
      args:
        type: object
        additionalProperties: false
        required: [series_matrix, output_dir]
        properties:
          series_matrix:
            type: string
            format: input-path
            description: GEO GSE*_series_matrix.txt or .txt.gz file.
          output_dir:
            type: string
            format: project-path
            description: Output directory for samplesheet.tsv and expression_matrix.tsv.
          platform_table:
            type: string
            format: input-path
            description: Optional GPL platform annotation TSV/CSV for probe-to-gene mapping.
          probe_col:
            type: string
            description: Optional probe ID column in platform_table (default ID or first column).
          gene_col:
            type: string
            description: Optional gene symbol column in platform_table.
          multi_gene_probes:
            type: string
            enum: [primary, drop]
            description: Handling of 'GENE1 /// GENE2' multi-mapped probes (default primary).
      approval: write
      output: ./schemas/geo-import-result.json
      timeoutSeconds: 300
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
    - name: coexpression
      description: >-
        Run WGCNA-style weighted gene co-expression network analysis on a bulk
        cohort (N >= 6, recommended >= 12): evaluate scale-free topology R^2 to
        select soft-threshold power beta, build signed/unsigned TOM hierarchical
        modules, merge correlated Module Eigengenes, rank intramodular kME hub
        genes, and test Module Eigengene-phenotype associations with BH FDR.
      run: [python, ./scripts/rna_coexpression.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Gene x sample expression matrix TSV/CSV (raw counts or normalized).
          output_dir:
            type: string
            format: project-path
            description: Directory for scale_free_fit.tsv, module_eigengenes.tsv, gene_modules_hubs.tsv, and summary plot.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV for Module Eigengene-trait correlations.
          sample_col:
            type: string
            description: Optional sample ID column in metadata.
          trait_cols:
            type: array
            items:
              type: string
            description: Optional clinical/phenotype trait columns in metadata to correlate with Module Eigengenes.
          top_genes:
            type: integer
            minimum: 50
            description: Number of top variable genes (by MAD) included in TOM construction (default 1500).
          network_type:
            type: string
            enum: [signed, unsigned]
            description: Co-expression similarity type (default signed).
          power:
            type: integer
            minimum: 0
            description: Soft-threshold power beta (0 = automatic scale-free R^2 selection, default 0).
          sft_r2_cutoff:
            type: number
            description: Scale-free topology R^2 threshold for automatic power selection (default 0.80).
          min_module_size:
            type: integer
            minimum: 5
            description: Minimum genes per co-expression module (default 25).
          merge_cut_height:
            type: number
            description: Merge modules with eigengene correlation > 1 - merge_cut_height (default 0.20).
      approval: write
      output: ./schemas/coexpression-result.json
      timeoutSeconds: 900
    - name: tme_deconv
      description: >-
        Estimate tumor microenvironment stromal/immune infiltration and purity
        (ESTIMATE-style StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity)
        plus non-negative least squares (NNLS) relative immune cell fractions
        across 10 canonical leukocyte lineages (or a custom signature matrix),
        with optional group differential infiltration testing.
      run: [python, ./scripts/rna_tme_deconv.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Gene x sample expression matrix TSV/CSV (HGNC/MGI gene symbols; map Ensembl IDs via rna_id_map first).
          output_dir:
            type: string
            format: project-path
            description: Directory for tme_cell_fractions.tsv, tme_microenvironment_scores.tsv, and summary plot.
          signature_matrix:
            type: string
            format: input-path
            description: Optional custom linear reference signature matrix TSV/CSV (genes x cell_types).
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV for group infiltration comparison.
          sample_col:
            type: string
            description: Optional sample ID column in metadata.
          group_col:
            type: string
            description: Optional condition/group column in metadata for Mann-Whitney/Kruskal + BH FDR tests.
      approval: write
      output: ./schemas/tme-deconv-result.json
      timeoutSeconds: 600
    - name: regulon
      description: >-
        Infer signed transcription factor and signaling pathway regulon
        activities from an expression matrix or DE result table using
        decoupler-style Univariate (ULM) or ridge-stabilized Multivariate (MLM)
        linear models with a built-in 14-pathway human/mouse footprint network
        or custom CollecTRI/DoRothEA/PROGENy regulons.
      run: [python, ./scripts/rna_regulon.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Expression matrix (genes x samples) or DE result table TSV/CSV.
          output_dir:
            type: string
            format: project-path
            description: Directory for regulon_activities.tsv, regulon_pvalues.tsv, regulon_fdr.tsv, and heatmap.
          network:
            type: string
            format: input-path
            description: Optional signed regulon TSV/CSV (source, target, optional mor/weight); uses built-in 14-pathway footprint resource if omitted.
          method:
            type: string
            enum: [ulm, mlm]
            description: Linear regression model (ulm or mlm, default ulm).
          min_targets:
            type: integer
            minimum: 3
            description: Minimum matched targets required per regulon (default 5).
          stat_col:
            type: string
            description: Optional statistic/log2FC column name when input is a DE table.
      approval: write
      output: ./schemas/regulon-result.json
      timeoutSeconds: 600
    - name: splicing
      description: >-
        Harmonize rMATS (SE, A5SS, A3SS, MXE, RI), SUPPA2 (.dpsi), and generic
        differential alternative splicing tables; compute exon lengths and flag
        frameshift_candidate (len % 3 != 0) vs in_frame_modular events; and
        summarize event classes and multi-hit spliced genes.
      run: [python, ./scripts/rna_splicing.py]
      args:
        type: object
        additionalProperties: false
        required: [splicing_inputs, output_dir]
        properties:
          splicing_inputs:
            type: array
            items:
              type: string
            description: One or more rMATS, SUPPA2, or generic differential splicing TSV/CSV tables.
          output_dir:
            type: string
            format: project-path
            description: Directory for splicing_significant_events.tsv, splicing_event_class_summary.tsv, multi_hit_genes, and volcano plot.
          dpsi_cutoff:
            type: number
            description: Minimum |dPSI| / |IncLevelDifference| threshold (default 0.10).
          fdr_cutoff:
            type: number
            description: Maximum FDR significance cutoff (default 0.05).
      approval: write
      output: ./schemas/splicing-result.json
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
expression + GMT ──────rna_pathway_score──▶ sample-level pathway activity matrix + heatmap
expression / DE ───────rna_regulon────────▶ ULM/MLM TF & pathway regulon activity t-scores + FDR
rMATS / SUPPA2 ────────rna_splicing───────▶ harmonized splicing events, frame impact & multi-hit genes
expression (+ traits) ─rna_coexpression───▶ WGCNA-style TOM modules, Module Eigengenes & kME hubs
symbol expression ─────rna_tme_deconv─────▶ ESTIMATE Stromal/Immune/Purity + NNLS immune cell fractions
```

## Hard red lines

1. **Never pass TPM/FPKM/VST to DESeq2 or edgeR**: verify `inferred_kind == "raw_counts"` via `rna_inspect` first.
2. **Check batch-condition collinearity**: if `rna_inspect` flags `confounded: true`, do not silently fit `~ batch + condition` (rank-deficient).
3. **Explicit background universe in ORA**: pass the full DE table (or explicit `--background` of expressed genes), never an unstated genome-wide default.
4. **Unfiltered ranked table in GSEA**: never pre-filter the DE table to `padj < 0.05` before calling `rna_gsea`.
5. **Replicate floor**: warn explicitly whenever any group has `n < 3` replicates.
