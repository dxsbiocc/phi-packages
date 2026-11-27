---
name: gwas-workflow
description: >-
  Statistical genetics and post-GWAS interpretation workflow in Python:
  standardize PLINK2/REGENIE/SAIGE/BOLT-LMM/GWAS Catalog summary statistics
  (lambda_GC, lambda_1000, lead loci, Manhattan & QQ plots); estimate SNP
  heritability h2, LDSC intercept, attenuation ratio, and cross-trait genetic
  correlation rg via block-jackknife LD Score Regression; run multi-causal
  SuSiE-RSS fine-mapping (or Wakefield ABF fallback) for per-SNP PIPs and
  purity-filtered 95% credible sets; align alleles and run FE-IVW / RE-DL
  cross-cohort GWAS meta-analysis with Cochran's Q, I^2 & EAF QC; and evaluate
  polygenic risk scores (PRS) for incremental delta-R2/AUC/C-index and quantile
  risk stratification. Use for GWAS, LDSC, heritability, genetic correlation,
  fine-mapping, SuSiE, credible set, meta-analysis, PRS, polygenic risk,
  Manhattan, 遗传力, 精细定位, 多基因风险评分, 荟萃分析.
phi:
  environment: plugin:gwas
  attachTo:
    - StatisticalGenetics
  scripts:
    - name: sumstats_harmonize
      description: >-
        Standardize GWAS summary statistics (PLINK2, REGENIE, SAIGE, BOLT-LMM,
        or GWAS Catalog), convert OR/LOG10P, enforce IUPAC alleles & SE > 0,
        filter MAF/INFO, compute lambda_GC & lambda_1000, extract clumped lead
        loci, and render a 2-panel Manhattan + QQ plot.
      run: [python, ./scripts/gwas_sumstats_harmonize.py]
      args:
        type: object
        additionalProperties: false
        required: [sumstats, output_dir]
        properties:
          sumstats:
            type: string
            format: input-path
            description: GWAS summary-statistics TSV/CSV (.gz supported).
          output_dir:
            type: string
            format: project-path
            description: Output directory for harmonized_sumstats.tsv, lead_loci.tsv, and gwas_manhattan_qq.png.
          n_cases:
            type: integer
            minimum: 1
            description: Optional case count for effective N (N_EFF) and case-control lambda_1000 calculation.
          n_controls:
            type: integer
            minimum: 1
            description: Optional control count for effective N (N_EFF) and case-control lambda_1000 calculation.
          min_maf:
            type: number
            minimum: 0
            description: Minimum minor allele frequency filter when EAF is present (default 0.01).
          min_info:
            type: number
            minimum: 0
            description: Minimum imputation INFO score filter when INFO is present (default 0.80).
          lead_window_kb:
            type: integer
            minimum: 1
            description: Minimum physical distance (kb) between independent lead loci (default 500).
          p_threshold:
            type: number
            description: Genome-wide significance P-value threshold for lead loci (default 5e-8).
      approval: write
      output: ./schemas/sumstats-harmonize-result.json
      timeoutSeconds: 600
    - name: ldsc_interpret
      description: >-
        Run Bulik-Sullivan (2015) heteroscedasticity-weighted LD Score
        Regression with genomic block jackknife SEs to estimate observed and
        liability-scale h2, LDSC intercept, attenuation ratio, and optional
        bivariate cross-trait genetic correlation rg.
      run: [python, ./scripts/gwas_ldsc_interpret.py]
      args:
        type: object
        additionalProperties: false
        required: [trait1_sumstats, output_dir]
        properties:
          trait1_sumstats:
            type: string
            format: input-path
            description: Primary trait GWAS summary-statistics TSV/CSV with SNP, Z (or BETA+SE), N, and optional l2.
          output_dir:
            type: string
            format: project-path
            description: Output directory for ldsc_summary.tsv and ldsc_binned_chi2_fit.png.
          ld_scores:
            type: string
            format: input-path
            description: Optional external LD score table containing SNP and l2 columns.
          trait2_sumstats:
            type: string
            format: input-path
            description: Optional second trait summary-statistics table for bivariate genetic correlation rg.
          sample_prev:
            type: number
            description: Optional sample case proportion (0..1) for Lee et al. (2011) liability-scale h2.
          pop_prev:
            type: number
            description: Optional population prevalence (0..1) for Lee et al. (2011) liability-scale h2.
          n_blocks:
            type: integer
            minimum: 10
            description: Number of contiguous genomic block-jackknife blocks (default 50, minimum 10).
      approval: write
      output: ./schemas/ldsc-interpret-result.json
      timeoutSeconds: 600
    - name: finemap
      description: >-
        Run SuSiE-RSS multi-causal statistical fine-mapping on single-locus
        summary statistics and regularized LD matrix R (or single-causal L=1
        Wakefield ABF fallback when ld_matrix is omitted), outputting per-SNP
        PIPs, purity-filtered 95% credible sets, and a locus PIP plot.
      run: [python, ./scripts/gwas_finemap.py]
      args:
        type: object
        additionalProperties: false
        required: [locus_sumstats, output_dir]
        properties:
          locus_sumstats:
            type: string
            format: input-path
            description: Single-locus summary-statistics TSV/CSV with SNP, POS, BETA+SE or Z, and N.
          output_dir:
            type: string
            format: project-path
            description: Output directory for finemap_pip.tsv, credible_sets.tsv, and finemap_pip_locus.png.
          ld_matrix:
            type: string
            format: input-path
            description: Optional M x M square LD correlation matrix (whitespace/TSV/CSV) matching locus_sumstats order.
          max_causal:
            type: integer
            minimum: 1
            maximum: 10
            description: Maximum number of causal variants L in SuSiE-RSS (default 5, range 1..10).
          coverage:
            type: number
            description: Target credible-set posterior probability coverage (default 0.95).
          min_abs_corr:
            type: number
            description: Minimum pairwise |r| purity filter for SuSiE credible sets (default 0.50).
      approval: write
      output: ./schemas/finemap-result.json
      timeoutSeconds: 600
    - name: meta_qc
      description: >-
        Align alleles across >=2 GWAS cohorts, flag cross-cohort EAF outliers,
        run Fixed-Effect IVW and DerSimonian-Laird Random-Effects meta-analysis
        with Cochran's Q, I^2, and direction strings, and render heterogeneity
        and EAF concordance diagnostics.
      run: [python, ./scripts/gwas_meta_qc.py]
      args:
        type: object
        additionalProperties: false
        required: [cohorts, output_dir]
        properties:
          cohorts:
            type: array
            minItems: 2
            items:
              type: string
              format: input-path
            description: Two or more cohort GWAS summary-statistics TSV/CSV files.
          output_dir:
            type: string
            format: project-path
            description: Output directory for gwas_meta_sumstats.tsv, gwas_meta_top_hits.tsv, and gwas_meta_heterogeneity.png.
          cohort_names:
            type: array
            items:
              type: string
            description: Optional labels corresponding to each cohort file.
          eaf_diff_max:
            type: number
            minimum: 0
            description: Maximum cross-cohort EAF discrepancy before flagging an EAF outlier (default 0.20).
      approval: write
      output: ./schemas/meta-qc-result.json
      timeoutSeconds: 600
    - name: prs_eval
      description: >-
        Standardize a polygenic risk score (prs_z), compare Clinical Covariates
        Only vs Covariates + PRS models (delta-R2 / Nagelkerke / liability R2,
        LRT P, delta-AUC or Cox PH delta-C-index), and compute OR/HR risk
        stratification across PRS quantiles vs Middle 40-60%.
      run: [python, ./scripts/gwas_prs_eval.py]
      args:
        type: object
        additionalProperties: false
        required: [cohort, prs_col, outcome_col, output_dir]
        properties:
          cohort:
            type: string
            format: input-path
            description: Individual-level TSV/CSV containing PRS, outcome, covariates, and optional survival time.
          prs_col:
            type: string
            description: Column name for the raw polygenic risk score.
          outcome_col:
            type: string
            description: Column name for the phenotype (continuous, binary 0/1, or event indicator).
          output_dir:
            type: string
            format: project-path
            description: Output directory for prs_quantile_risk.tsv, prs_model_comparison.tsv, and prs_stratification_plot.png.
          covariates:
            type: array
            items:
              type: string
            description: Optional covariate column names (e.g. age, sex, PC1, PC2).
          time_col:
            type: string
            description: Optional follow-up time column for Cox proportional hazards HR and C-index evaluation.
          pop_prev:
            type: number
            description: Optional population prevalence (0..1) for Lee et al. (2012) liability-scale R2.
      approval: write
      output: ./schemas/prs-eval-result.json
      timeoutSeconds: 600
---

# Statistical Genetics & Post-GWAS Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for
harmonization conventions, LDSC intercept/liability formulas, SuSiE-RSS IBSS
regularization and purity criteria, meta-analysis heterogeneity QC, and PRS
incremental evaluation guidelines.

## Chain

```text
Raw Genotypes / Summary Stats
  ├── (Upstream nf-core wrappers: plink2, regenie, gcta, beagle5, shapeit5)
  ├── gwas_sumstats_harmonize ─▶ standardize BETA/SE/alleles/LOG10P, λ_GC & λ_1000, lead loci
  ├── gwas_ldsc_interpret ─────▶ block-jackknife LDSC h2 (obs & liab), intercept, cross-trait rg
  ├── gwas_finemap ────────────▶ SuSiE-RSS (or ABF L=1) PIPs & purity-filtered 95% credible sets
  ├── gwas_meta_qc ────────────▶ cross-cohort allele/EAF alignment, FE-IVW & RE-DL, Q & I^2
  └── gwas_prs_eval ───────────▶ standardized PRS ΔR2/ΔAUC/ΔC-index & quantile OR/HR stratification
```
