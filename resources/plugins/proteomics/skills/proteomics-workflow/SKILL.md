---
name: proteomics-workflow
description: >-
  Downstream quantitative MS-based proteomics and PTM/phosphoproteomics
  analysis in Python: inspect DIA-NN / MaxQuant / FragPipe / generic protein
  and PTM intensity matrices (catching contaminants, decoys, raw-vs-log2 scale,
  missingness, and batch-condition collinearity); filter contaminants, log2
  transform, normalize (median/quantile), impute missing values (Perseus-style
  MinProb MNAR or kNN MAR), and run PCA/CV QC; run empirical-Bayes moderated
  t-test / Welch / Mann-Whitney differential abundance; normalize PTM/phospho
  sites against total protein abundance and run Kinase-Substrate Enrichment
  Analysis (KSEA); run offline detected-proteome ORA and pre-ranked GSEA
  against GMT gene sets; and build co-abundance correlation modules and hub
  rankings. Use for proteomics, 蛋白质组, 蛋白组, 磷酸化, phosphoproteomics,
  PTM, DIA-NN, MaxQuant, FragPipe, LFQ, TMT, KSEA, or protein differential
  abundance. Raw spectrum search/quantification runs upstream (quantms /
  DIA-NN / MaxQuant / FragPipe).
phi:
  environment: plugin:proteomics
  attachTo:
    - Proteomics
  scripts:
    - name: inspect
      description: >-
        Inspect a protein-group or PTM/phosphosite intensity matrix (MaxQuant
        proteinGroups, DIA-NN pg_matrix, FragPipe combined_protein, or generic
        TSV/CSV) and optional sample metadata: auto-detect format, contaminant/
        reverse counts, raw vs log2 scale, per-sample missingness, sample
        alignment, replicate counts, and batch-condition collinearity. Call
        first before QC or differential analysis.
      run: [python, ./scripts/prot_inspect.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: Protein-group or PTM intensity matrix TSV/CSV.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV.
          id_col:
            type: string
            description: Optional protein/site identifier column (auto-detected if omitted).
          sample_col:
            type: string
            description: Metadata column holding sample IDs (default first column).
          group_col:
            type: string
            description: Primary condition column in metadata (default condition).
          batch_col:
            type: string
            description: Optional batch column to test for collinearity with group_col.
      approval: read
      timeoutSeconds: 300
    - name: rollup
      description: >-
        Aggregate peptide or precursor intensity tables (DIA-NN pr_matrix,
        MaxQuant peptides.txt, FragPipe combined_peptide, or QUANTIX/generic
        peptide tables) to log2 protein-level abundances using reference-scaled
        rollup (rrollup), Top-3 (Hi3), median, or linear sum, with unique-only
        or Occam's razor shared-peptide assignment and contaminant removal.
      run: [python, ./scripts/prot_rollup.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Peptide or precursor intensity matrix TSV/CSV.
          peptide_col:
            type: string
            description: Optional peptide/precursor identifier column (auto-detected if omitted).
          protein_col:
            type: string
            description: Optional protein/gene identifier column (auto-detected if omitted).
          shared_peptides:
            type: string
            enum: [unique_only, razor]
            description: Shared-peptide assignment policy (default unique_only; razor assigns to most-supported protein).
          method:
            type: string
            enum: [rrollup, top3, median, sum]
            description: Protein rollup strategy (default rrollup reference-scaled log2 rollup).
          min_peptides:
            type: integer
            minimum: 1
            description: Minimum peptides required to quantify a protein (default 2).
          output_dir:
            type: string
            format: project-path
            description: Directory for rolled-up log2 protein matrix, peptide count table, and histogram.
      approval: write
      timeoutSeconds: 600
    - name: qc_impute
      description: >-
        Remove contaminants/decoys, log2-transform raw intensities (if needed),
        filter proteins by valid-value completeness (overall or within at least
        one group), apply median or quantile normalization, impute missing
        values (Perseus-style left-censored MinProb MNAR, kNN MAR, or none),
        and generate processed matrix plus missingness, intensity-distribution,
        and PCA QC diagnostics.
      run: [python, ./scripts/prot_qc_impute.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Protein or PTM intensity matrix TSV/CSV.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV (for group-aware filtering and PCA coloring).
          id_col:
            type: string
            description: Optional protein/site identifier column (auto-detected if omitted).
          sample_col:
            type: string
            description: Metadata sample ID column (default first column).
          group_col:
            type: string
            description: Metadata condition column (default condition).
          min_valid_frac:
            type: number
            description: Minimum fraction of non-missing values required overall or in >=1 group (default 0.5).
          filter_mode:
            type: string
            enum: [in_at_least_one_group, overall]
            description: Apply min_valid_frac per group (default when metadata present) or globally.
          normalization:
            type: string
            enum: [median, quantile, none]
            description: Sample normalization method on log2 scale (default median).
          imputation:
            type: string
            enum: [minprob, knn, none]
            description: Missing-value imputation method (default minprob for MNAR; knn for MAR; none to skip).
          downshift:
            type: number
            description: MinProb mean downshift in SD units (default 1.8).
          width:
            type: number
            description: MinProb SD width scaling factor (default 0.3).
          knn_k:
            type: integer
            minimum: 1
            description: Number of neighbours for kNN imputation (default 5).
          seed:
            type: integer
            description: Random seed for reproducible MinProb imputation (default 42).
          output_dir:
            type: string
            format: project-path
            description: Directory for processed matrix, sample QC table, and QC figures.
      approval: write
      timeoutSeconds: 600
    - name: de
      description: >-
        Differential protein abundance between two conditions using
        empirical-Bayes moderated t-test (limma-style variance shrinkage),
        Welch t-test, or Mann-Whitney U with Benjamini-Hochberg FDR, writing
        the full DE table plus volcano and MA plots. Expects log2-scale
        normalized intensities (output of prot_qc_impute).
      run: [python, ./scripts/prot_de.py]
      args:
        type: object
        additionalProperties: false
        required: [input, metadata, case, control, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Normalized log2 protein intensity matrix TSV/CSV (from prot_qc_impute).
          metadata:
            type: string
            format: input-path
            description: Sample metadata TSV/CSV.
          id_col:
            type: string
            description: Optional protein identifier column (default first column).
          sample_col:
            type: string
            description: Metadata sample ID column (default first column).
          group_col:
            type: string
            description: Metadata condition column (default condition).
          case:
            type: string
            description: Case / treatment group label in group_col.
          control:
            type: string
            description: Control / reference group label in group_col.
          method:
            type: string
            enum: [moderated_t, welch, mannwhitney]
            description: Statistical test (default moderated_t).
          fdr:
            type: number
            description: BH FDR cutoff for calling regulated proteins (default 0.05).
          log2fc:
            type: number
            description: Minimum absolute log2 fold-change cutoff (default 1.0).
          output_dir:
            type: string
            format: project-path
            description: Directory for DE table, volcano plot, and MA plot.
      approval: write
      timeoutSeconds: 600
    - name: ptm
      description: >-
        PTM / phosphoproteomics site-level differential analysis with optional
        total-protein abundance normalization (site_log2 - protein_log2) and
        Kinase-Substrate Enrichment Analysis (KSEA z-score + BH FDR) against a
        kinase-substrate mapping table.
      run: [python, ./scripts/prot_ptm.py]
      args:
        type: object
        additionalProperties: false
        required: [ptm_matrix, metadata, case, control, output_dir]
        properties:
          ptm_matrix:
            type: string
            format: input-path
            description: Log2 PTM/phosphosite intensity matrix TSV/CSV (with site_id and optional protein/gene column).
          protein_matrix:
            type: string
            format: input-path
            description: Optional matched log2 total-protein matrix for protein-abundance normalization.
          metadata:
            type: string
            format: input-path
            description: Sample metadata TSV/CSV.
          site_col:
            type: string
            description: Column holding PTM site IDs (default first column, e.g. site_id).
          protein_col:
            type: string
            description: Optional column in ptm_matrix linking each site to its parent protein/gene ID.
          sample_col:
            type: string
            description: Metadata sample ID column (default first column).
          group_col:
            type: string
            description: Metadata condition column (default condition).
          case:
            type: string
            description: Case / treatment condition label.
          control:
            type: string
            description: Control / reference condition label.
          kinase_table:
            type: string
            format: input-path
            description: Optional TSV/CSV with kinase and substrate (or site_id) columns for KSEA.
          min_substrates:
            type: integer
            minimum: 2
            description: Minimum quantified substrates required to score a kinase (default 3).
          fdr:
            type: number
            description: BH FDR cutoff for regulated sites and kinases (default 0.05).
          log2fc:
            type: number
            description: Minimum absolute log2FC cutoff for regulated PTM sites (default 1.0).
          output_dir:
            type: string
            format: project-path
            description: Directory for PTM DE table, KSEA table, and figures.
      approval: write
      timeoutSeconds: 600
    - name: enrichment
      description: >-
        Offline pathway enrichment from a full proteomics DE table against a
        GMT gene-set file: hypergeometric ORA strictly against the detected-
        proteome universe (never a full-genome default) and/or deterministic
        pre-ranked GSEA (ES/NES/FDR + leading-edge proteins) plus dotplot.
      run: [python, ./scripts/prot_enrichment.py]
      args:
        type: object
        additionalProperties: false
        required: [de_table, gmt, output_dir]
        properties:
          de_table:
            type: string
            format: input-path
            description: Full unfiltered DE table from prot_de.
          gmt:
            type: string
            format: input-path
            description: GMT gene-set file (gene symbols matching de_table IDs).
          mode:
            type: string
            enum: [both, ora, gsea]
            description: Which enrichment method(s) to run (default both).
          direction:
            type: string
            enum: [both, up, down]
            description: Foreground direction for ORA (default both).
          fdr:
            type: number
            description: DE and pathway BH FDR cutoff (default 0.05).
          log2fc:
            type: number
            description: DE |log2FC| cutoff for ORA foreground (default 1.0).
          min_set_size:
            type: integer
            minimum: 2
            description: Minimum overlap with detected proteome (default 5).
          max_set_size:
            type: integer
            minimum: 5
            description: Maximum overlap with detected proteome (default 500).
          permutations:
            type: integer
            minimum: 100
            description: GSEA permutation count (default 500).
          seed:
            type: integer
            description: Random seed for GSEA permutations (default 42).
          output_dir:
            type: string
            format: project-path
            description: Directory for ORA/GSEA tables and enrichment dotplot.
      approval: write
      timeoutSeconds: 1200
    - name: network
      description: >-
        Build a protein co-abundance correlation network across top variable or
        regulated proteins, extract connected modules and weighted hub scores,
        and write edge list, hub ranking, and clustered correlation heatmap.
      run: [python, ./scripts/prot_network.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Normalized log2 protein intensity matrix TSV/CSV.
          de_table:
            type: string
            format: input-path
            description: Optional prot_de output to prioritize regulated proteins.
          id_col:
            type: string
            description: Optional protein identifier column (default first column).
          top_n:
            type: integer
            minimum: 10
            description: Maximum number of top variable/regulated proteins to include (default 80).
          method:
            type: string
            enum: [spearman, pearson]
            description: Correlation metric (default spearman).
          min_abs_corr:
            type: number
            description: Minimum |r| threshold to retain an edge (default 0.7).
          fdr:
            type: number
            description: BH FDR cutoff for correlation edges (default 0.05).
          output_dir:
            type: string
            format: project-path
            description: Directory for edge table, hub table, and co-abundance heatmap.
      approval: write
      timeoutSeconds: 600
---

# Proteomics — quantitative MS & PTM downstream chain

Read [references/analysis-playbook.md](references/analysis-playbook.md) before
choosing normalization, imputation (`minprob` vs `knn`), and PTM normalization
strategies.

## Chain

```
Raw MS (.raw/.d/.mzML) ──(quantms / DIA-NN / MaxQuant / FragPipe)──▶ peptide / protein / PTM matrix
peptide matrix    ──prot_rollup─────▶ unique/razor peptide-to-protein rollup (rrollup/top3/median/sum)
matrix + metadata ──prot_inspect────▶ detect format, contaminants, scale, missingness & design
                  ──prot_qc_impute──▶ filter contaminants, log2, normalize, impute & QC/PCA
                  ──prot_de─────────▶ empirical-Bayes moderated t-test + volcano/MA
PTM + total prot  ──prot_ptm────────▶ protein-normalized PTM DE + KSEA kinase activity
full DE + GMT     ──prot_enrichment─▶ detected-proteome ORA + pre-ranked GSEA
processed matrix  ──prot_network────▶ co-abundance modules & hub proteins
```

## Hard red lines

1. **Remove contaminants & decoys first**: never keep `CON__`, `REV__`, `Potential contaminant == "+"`, or `Reverse == "+"` rows in downstream statistics.
2. **Never run t-tests on raw linear intensities**: verify `log2` scale via `prot_inspect` and `prot_qc_impute`.
3. **Match imputation to missingness**: `minprob` for MNAR low-abundance dropouts (LFQ/DIA), `knn` for MAR dropouts (TMT), `none` if nearly complete.
4. **Normalize PTM sites by total protein whenever available**: without `--protein_matrix`, state explicitly that site fold-changes confound abundance and stoichiometry.
5. **Detected-proteome universe in ORA**: always use the full quantified proteome from `prot_de`, never an unstated genome-wide background.
