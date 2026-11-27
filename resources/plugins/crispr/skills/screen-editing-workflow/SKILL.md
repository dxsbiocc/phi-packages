---
name: screen-editing-workflow
description: >-
  Pooled CRISPR functional & chemogenetic screen QC, alpha-RRA & DrugZ hit
  calling, and post-editing amplicon outcome quantification: audit count
  matrices from mageck count or crisprcleanr (read depth, zero-dropout %, Lorenz
  Gini index, DESeq median normalization, PCA/correlation, and essential-vs-nonessential
  NNMD & ROC-AUC); call positive and negative selection gene hits via alpha-RRA;
  call drug-vs-vehicle synthetic-lethal sensitizers and resistance suppressors
  via Collimore & Hart (2019) DrugZ empirical-Bayes normZ; and quantify NHEJ
  frameshift KO %, in-frame indels, clean vs imperfect HDR %, and CBE/ABE
  target-vs-bystander editing from amplicon allele tables/FASTAs or CRISPResso2
  outputs. Use for MAGeCK, DrugZ, synthetic lethality, CRISPRcleanR, CRISPResso2,
  Gini, NNMD, RRA, indel, frameshift, HDR, base editing, 筛选质控, 富集筛选, 合成致死, 编辑效率定量.
phi:
  environment: plugin:crispr
  attachTo:
    - Crispr
  scripts:
    - name: screen_qc
      description: >-
        Evaluate pooled CRISPR screen raw count matrix quality (compatible with
        mageck count .count.txt and generic TSV/CSV), reporting median reads per
        sgRNA, zero-dropout %, Lorenz Gini index, DESeq median-of-ratios size
        factors, sample Pearson correlation & PCA, and essential vs nonessential
        NNMD & ROC-AUC when control gene lists are provided.
      run: [python, ./scripts/crispr_screen_qc.py]
      args:
        type: object
        additionalProperties: false
        required: [counts, output_dir]
        properties:
          counts:
            type: string
            format: input-path
            description: sgRNA count matrix TSV/CSV (e.g. mageck count .count.txt or crisprcleanr normalized counts).
          output_dir:
            type: string
            format: project-path
            description: Directory for screen_sample_qc.tsv and screen_qc_dashboard.png.
          sgrna_col:
            type: string
            description: sgRNA ID column name (auto-detected if omitted).
          gene_col:
            type: string
            description: Target gene column name (auto-detected if omitted).
          metadata:
            type: string
            format: input-path
            description: Optional sample sheet TSV/CSV with sample, condition, and timepoint columns.
          essential_genes:
            type: string
            format: input-path
            description: Optional newline-delimited or TSV list of core essential genes (Hart/DepMap).
          nonessential_genes:
            type: string
            format: input-path
            description: Optional newline-delimited or TSV list of non-essential / NTC genes.
      approval: write
      output: ./schemas/screen-qc-result.json
      timeoutSeconds: 300
    - name: screen_hits
      description: >-
        Call positive (enrichment) and negative (depletion) selection gene hits
        from a pooled CRISPR sgRNA count matrix using count normalization,
        moderated sgRNA variance z-scores, MAGeCK-inspired alpha-RRA scoring,
        seeded sgRNA-multiplicity permutation p-values, and BH FDR.
      run: [python, ./scripts/crispr_screen_hits.py]
      args:
        type: object
        additionalProperties: false
        required: [counts, treatment_cols, control_cols, output_dir]
        properties:
          counts:
            type: string
            format: input-path
            description: sgRNA raw or CNV-normalized count matrix TSV/CSV.
          treatment_cols:
            type: array
            items:
              type: string
            description: Treatment/endpoint sample column names in counts.
          control_cols:
            type: array
            items:
              type: string
            description: Control/baseline/T0 sample column names in counts.
          output_dir:
            type: string
            format: project-path
            description: Directory for sgrna_summary.tsv, gene_hits_summary.tsv, and screen_hits_volcano_rank.png.
          sgrna_col:
            type: string
            description: sgRNA ID column name (auto-detected if omitted).
          gene_col:
            type: string
            description: Target gene column name (auto-detected if omitted).
          paired:
            type: boolean
            description: Compute paired replicate log2FC when treatment_cols and control_cols are matched 1:1 (default false).
          norm_method:
            type: string
            enum: [median, control_sgrna, total]
            description: Count normalization strategy (default median).
          control_gene_prefix:
            type: string
            description: Comma-separated prefixes for non-targeting control genes when norm_method=control_sgrna (default NTC,NonTargeting).
          alpha_cutoff:
            type: number
            description: Quantile cutoff alpha for alpha-RRA rank aggregation (default 0.10).
          n_permutations:
            type: integer
            minimum: 200
            description: Number of seeded empirical permutations per gene size for RRA p-value calibration (default 2000).
      approval: write
      output: ./schemas/screen-hits-result.json
      timeoutSeconds: 600
    - name: drugz
      description: >-
        Chemogenetic drug-vs-vehicle pooled CRISPR screen analysis implementing
        the Collimore & Hart (2019) DrugZ algorithm: 10M RPM normalization,
        abundance-ordered empirical-Bayes monotonic sliding-window variance
        estimation, paired/unpaired replicate z-scores, gene-level sumZ/normZ,
        and BH FDR for synthetic-lethal sensitizers and drug-resistance suppressors.
      run: [python, ./scripts/crispr_drugz.py]
      args:
        type: object
        additionalProperties: false
        required: [counts, drug_cols, control_cols, output_dir]
        properties:
          counts:
            type: string
            format: input-path
            description: sgRNA count matrix TSV/CSV (e.g. mageck count .count.txt).
          drug_cols:
            type: array
            items:
              type: string
            description: Drug-treated replicate column names in counts.
          control_cols:
            type: array
            items:
              type: string
            description: Vehicle/DMSO-treated control replicate column names in counts.
          output_dir:
            type: string
            format: project-path
            description: Directory for drugz_sgrna_zscores.tsv, drugz_gene_summary.tsv, and drugz_normz_waterfall.png.
          sgrna_col:
            type: string
            description: sgRNA ID column name (auto-detected if omitted).
          gene_col:
            type: string
            description: Target gene column name (auto-detected if omitted).
          unpaired:
            type: boolean
            description: Force unpaired group-mean DrugZ comparison even when drug_cols and control_cols have equal length (default false).
          pseudocount:
            type: number
            description: Stabilizing pseudocount added to 10M-normalized RPM reads before log2 fold-change (default 5.0).
          window_size:
            type: integer
            minimum: 25
            description: Abundance-ordered sliding window size for empirical-Bayes MAD variance estimation (default 500).
          remove_genes:
            type: array
            items:
              type: string
            description: Optional core-essential genes to exclude before variance estimation so extreme pan-essential dropout does not inflate background variance.
      approval: write
      output: ./schemas/drugz-result.json
      timeoutSeconds: 300
    - name: editing_quant
      description: >-
        Quantify post-editing genome modification outcomes from raw amplicon
        allele tables/FASTAs or nf-core/modules/crispresso2 allele frequency
        tables, distinguishing Unmodified WT, productive Frameshift Indel KO,
        In-frame Indel, Pure Substitution, clean HDR, and Imperfect HDR, plus
        indel size spectrum and protospacer CBE/ABE target vs bystander rates.
      run: [python, ./scripts/crispr_editing_quant.py]
      args:
        type: object
        additionalProperties: false
        required: [amplicon_seq, guide_seq, output_dir]
        properties:
          alleles:
            type: string
            format: input-path
            description: Mode A input — TSV/CSV with sequence + count columns (and optional sample), or amplicon FASTA.
          crispresso_table:
            type: string
            format: input-path
            description: Mode B input — Alleles_frequency_table_around_sgRNA_*.txt TSV produced by nf-core/modules/crispresso2.
          amplicon_seq:
            type: string
            description: Wild-type reference amplicon DNA sequence (5'->3').
          guide_seq:
            type: string
            description: 20-23 nt protospacer guide sequence located inside amplicon_seq (either strand).
          output_dir:
            type: string
            format: project-path
            description: Directory for editing_allele_summary.tsv, indel_size_distribution.tsv, protospacer_base_matrix.tsv, and editing_outcomes.png.
          hdr_seq:
            type: string
            description: Optional expected HDR donor amplicon sequence (5'->3') for knock-in vs imperfect-HDR classification.
          editor_mode:
            type: string
            enum: [nuclease, CBE, ABE]
            description: Editing modality determining base-conversion and product-purity metrics (default nuclease).
          window_bp:
            type: integer
            minimum: 2
            description: Half-window in bp around the cleavage site for editing classification (default 10).
      approval: write
      output: ./schemas/editing-quant-result.json
      timeoutSeconds: 300
---

# Pooled Screen QC, Hit Calling & Post-Editing Quantification Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for Gini/NNMD thresholds, CNV correction routing, $\alpha$-RRA equations, and frameshift/bystander quantification rules.

## Workflow Chain

```
Screen FASTQ + library ──nf-core/modules/mageck-count──────────▶ Raw sgRNA count matrix (.count.txt)
Cancer lethality screen──nf-core/modules/crisprcleanr-normalize▶ CNV-corrected sgRNA count matrix
Count matrix           ──crispr_screen_qc──────────────────────▶ Gini, zero-dropout, PCA, NNMD & ROC-AUC
Count matrix           ──crispr_screen_hits────────────────────▶ sgRNA log2FC + gene alpha-RRA + BH FDR
Drug vs DMSO counts    ──crispr_drugz──────────────────────────▶ Empirical-Bayes normZ synthetic lethality & resistance
Amplicon FASTQ         ──nf-core/modules/crispresso2───────────▶ Alleles_frequency_table_around_sgRNA_*.txt
Allele table / CRISPResso2 ──crispr_editing_quant──────────────▶ Frameshift KO %, HDR %, CBE/ABE bystander %
```

## Hard Red Lines

1. **Screen bottlenecking & library representation**: flag plasmid/T0 samples with Gini `> 0.15` (and any sample `> 0.25` / severe `> 0.30`) or median `< 100` reads/sgRNA (`< 300` caution). Verify essential vs non-essential separation (`NNMD < -1.25`) before trusting negative-selection hits.
2. **CNV false positives in cancer cell lines**: Cas9 cleavage at copy-number-amplified loci causes gene-independent DNA-damage lethality. Run `nf-core/modules/crisprcleanr-normalize` prior to `crispr_screen_hits` when profiling aneuploid cancer lines.
3. **Frameshift vs total indel distinction**: in KO validation, report productive **frameshift indel %** (`|ins - del| % 3 != 0`) separately from **in-frame indel %** (`|ins - del| % 3 == 0`).
4. **Base-editor product purity**: in `CBE`/`ABE` runs, always quantify **bystander edits** within protospacer positions `4..8` and **indel byproduct %** alongside target conversion %.
