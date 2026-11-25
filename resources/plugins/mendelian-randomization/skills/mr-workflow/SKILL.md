---
name: mr-workflow
description: >-
  Genetic epidemiology, GWAS summary-statistics QC, Two-Sample Mendelian
  Randomization (STROBE-MR), and Bayesian colocalization workflow in Python:
  standardize GWAS tables (OR-to-logOR, allele cleaning, genomic inflation
  lambda_GC, Manhattan & QQ plots, physical window clumping); harmonize exposure
  instruments against outcome GWAS (strand/allele flipping, palindromic EAF
  filtering, F>=10 weak-instrument & Steiger checks); run IVW (random/fixed),
  Wald ratio, MR-Egger + intercept pleiotropy, Weighted Median, Weighted Mode,
  Cochran's Q, outlier & Leave-One-Out diagnostics; and compute Wakefield ABF
  colocalization (PP.H0-PP.H4 + 95% credible set). Use for 孟德尔随机化,
  Mendelian randomization, TwoSampleMR, GWAS, Manhattan, 曼哈顿图, IVW,
  MR-Egger, 共定位, coloc, cis-pQTL, or cis-eQTL causal inference.
phi:
  environment: plugin:mr
  attachTo:
    - MendelianRandomization
  scripts:
    - name: gwas_qc
      description: >-
        Standardize a GWAS/QTL summary-statistics table (TSV/CSV/.gz; auto-maps
        SNP, CHR, POS, EA, OA, BETA/OR, SE, P, EAF, N; converts OR to log(OR)),
        compute genomic inflation factor lambda_GC, extract distance-clumped
        lead loci, and draw a 2-panel Manhattan + QQ plot.
      run: [python, ./scripts/mr_gwas_qc.py]
      args:
        type: object
        additionalProperties: false
        required: [gwas, output_dir]
        properties:
          gwas:
            type: string
            format: input-path
            description: GWAS or QTL summary-statistics TSV/CSV (.gz supported).
          p_cutoff:
            type: number
            description: P-value threshold for lead-variant extraction (default 5e-8).
          clump_kb:
            type: integer
            minimum: 50
            description: Physical clumping window in kb around lead variants (default 10000 for genome-wide IVs; use 500 for locus listing).
          min_maf:
            type: number
            description: Optional minimum minor allele frequency filter when EAF is present (default 0.01).
          sample_size:
            type: integer
            description: Optional cohort sample size N if absent from the GWAS columns.
          output_dir:
            type: string
            format: project-path
            description: Output directory for standardized GWAS, clumped leads, and Manhattan/QQ figure.
      approval: write
      output: ./schemas/gwas-qc-result.json
      timeoutSeconds: 600
    - name: harmonize
      description: >-
        Harmonize exposure genetic instruments against an outcome GWAS summary
        table: select/clump exposure SNPs below p_cutoff, enforce F >= min_f,
        align effect/other alleles (including unambiguous strand flips), filter
        ambiguous palindromic A/T & C/G SNPs by EAF margin, and optionally apply
        Steiger directionality filtering.
      run: [python, ./scripts/mr_harmonize.py]
      args:
        type: object
        additionalProperties: false
        required: [exposure, outcome, output_dir]
        properties:
          exposure:
            type: string
            format: input-path
            description: Standardized or raw exposure GWAS/QTL summary table (or pre-clumped IV table).
          outcome:
            type: string
            format: input-path
            description: Standardized or raw outcome GWAS summary table.
          p_cutoff:
            type: number
            description: Exposure instrument P-value cutoff (default 5e-8; relax to 1e-5 only when explicitly justified).
          clump_kb:
            type: integer
            minimum: 0
            description: Physical window clumping in kb on exposure SNPs when CHR/POS are present (default 10000; 0 disables if already LD-clumped).
          min_f:
            type: number
            description: Minimum per-SNP F-statistic beta^2/se^2 (default 10.0).
          palindromic_eaf_margin:
            type: number
            description: Drop palindromic SNPs with |EAF - 0.5| < margin (default 0.08, i.e. drop 0.42..0.58).
          steiger_filter:
            type: boolean
            description: Drop individual SNPs explaining more variance in outcome than exposure when N is available.
          output_dir:
            type: string
            format: project-path
            description: Output directory for harmonized_instruments.tsv.
      approval: write
      output: ./schemas/harmonize-result.json
      timeoutSeconds: 300
    - name: estimate
      description: >-
        Run Two-Sample Mendelian Randomization on harmonized_instruments.tsv:
        multiplicative random-effects IVW, fixed-effects IVW (or Wald ratio for
        1 SNP), MR-Egger + intercept pleiotropy test, Weighted Median, Simple &
        Weighted Mode, Cochran's Q & I^2 heterogeneity, residual outlier flags,
        Leave-One-Out sensitivity, global Steiger test, and 4-panel MR plots.
      run: [python, ./scripts/mr_estimate.py]
      args:
        type: object
        additionalProperties: false
        required: [harmonized, output_dir]
        properties:
          harmonized:
            type: string
            format: input-path
            description: Harmonized exposure-outcome table from mr_harmonize.
          exposure_name:
            type: string
            description: Exposure label for tables and plots (default Exposure).
          outcome_name:
            type: string
            description: Outcome label for tables and plots (default Outcome).
          outcome_binary:
            type: boolean
            description: Report exponentiated Odds Ratios (OR, 95% CI) alongside log-odds beta for binary disease outcomes.
          n_boot:
            type: integer
            minimum: 200
            description: Parametric bootstrap draws for Weighted Median and Mode SEs (default 1000).
          seed:
            type: integer
            description: Random seed for bootstrap reproducibility (default 42).
          output_dir:
            type: string
            format: project-path
            description: Output directory for mr_results.tsv, mr_singlesnp_loo.tsv, and mr-diagnostics.png.
      approval: write
      output: ./schemas/estimate-result.json
      timeoutSeconds: 300
    - name: coloc
      description: >-
        Run exact single-causal-variant Bayesian colocalization (Wakefield 2009
        ABF / Giambartolomei coloc.abf) across a shared genomic region between
        two GWAS/QTL summary tables: outputs posterior probabilities PP.H0-PP.H4,
        per-SNP conditional SNP.PP.H4, the 95% credible set, and a regional
        mirror plot.
      run: [python, ./scripts/mr_coloc.py]
      args:
        type: object
        additionalProperties: false
        required: [trait1, trait2, output_dir]
        properties:
          trait1:
            type: string
            format: input-path
            description: First locus summary-statistics table (e.g. cis-eQTL/pQTL or exposure GWAS).
          trait2:
            type: string
            format: input-path
            description: Second locus summary-statistics table (e.g. disease outcome GWAS).
          chrom:
            type: string
            description: Optional chromosome filter (e.g. '1' or 'chr1').
          start:
            type: integer
            description: Optional region start coordinate (bp).
          end:
            type: integer
            description: Optional region end coordinate (bp).
          trait1_type:
            type: string
            enum: [quant, cc]
            description: Trait 1 scale — 'quant' (prior W=0.15^2) or 'cc' case-control log-OR (prior W=0.20^2; default quant).
          trait2_type:
            type: string
            enum: [quant, cc]
            description: Trait 2 scale — 'quant' or 'cc' (default cc).
          p1:
            type: number
            description: Prior probability a SNP is associated with trait 1 only (default 1e-4).
          p2:
            type: number
            description: Prior probability a SNP is associated with trait 2 only (default 1e-4).
          p12:
            type: number
            description: Prior probability a SNP is causal for both traits (default 1e-5).
          output_dir:
            type: string
            format: project-path
            description: Output directory for coloc_snps.tsv and coloc-locus.png.
      approval: write
      output: ./schemas/coloc-result.json
      timeoutSeconds: 300
---

# GWAS, Colocalization & Mendelian Randomization Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for
STROBE-MR reporting rules, palindromic SNP handling, pleiotropy diagnostics, and
colocalization priors.

## Chain

```text
Exposure GWAS + Outcome GWAS
  ├── mr_gwas_qc ───▶ standardize BETA/SE/alleles, check λ_GC, Manhattan & QQ plot
  ├── mr_harmonize ─▶ clump/select IVs (F ≥ 10), align alleles, drop ambiguous palindromes
  ├── mr_estimate ──▶ IVW + MR-Egger + Weighted Median/Mode + Cochran Q + LOO + 4-panel plots
  └── mr_coloc ─────▶ cis-locus Wakefield ABF Bayesian colocalization (PP.H0..PP.H4 + 95% CS)
```
