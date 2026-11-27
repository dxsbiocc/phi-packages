---
name: pgx-workflow
description: >-
  Pharmacogenomics, high-throughput drug screening, cross-cohort sensitivity
  imputation, and combination synergy analysis: score HTS plate QC (Zhang's
  Z'-factor, SSMD, edge effects, Tukey B-score, MAD hit calling); fit
  4-parameter log-logistic Hill dose-response curves (EC50, absolute IC50, E0,
  Einf, Hill slope, normalized AUC/AAC, R^2, and beyond-tested-range censoring
  flags); score 2D drug-combination matrices with ZIP, Bliss, Loewe, and HSA
  synergy models; discover mutation/CNV or gene-expression biomarkers; and
  translate cell-line drug response models to clinical cohorts via RidgeCV. Use
  for pharmacogenomics, HTS plate QC, Z'-factor, B-score, IC50, EC50, AUC,
  drug screening, oncoPredict, 药敏分析, 剂量反应曲线, 药物联用协同,
  ZIP/Bliss/Loewe, GDSC/DepMap 药敏标志物, 患者药敏预测.
phi:
  environment: plugin:pgx
  attachTo:
    - Pharmacogenomics
  scripts:
    - name: dose_response
      description: >-
        Fit 4-parameter log-logistic Hill dose-response curves across compounds
        and cell lines/samples from a viability table, reporting relative EC50,
        absolute IC50 (y=0.5), E0, Einf (Emax viability), Hill slope, normalized
        trapezoidal AUC & AAC, R^2, and beyond_tested_range censoring flags.
      run: [python, ./scripts/pgx_dose_response.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Long-format dose-viability TSV/CSV with dose and viability columns (plus optional compound and sample columns).
          output_dir:
            type: string
            format: project-path
            description: Directory for fitted parameter table and dose-response curve plot.
          compound_col:
            type: string
            description: Compound/drug column name (default drug or compound if present).
          sample_col:
            type: string
            description: Cell line / sample column name (default sample or cell_line if present).
          dose_col:
            type: string
            description: Concentration/dose column name (default dose or conc, must be > 0; zero-dose vehicle rows are used for E0 normalization if present).
          viability_col:
            type: string
            description: Viability column name (default viability; auto-detects 0..1 fraction vs 0..100 percent).
          viability_scale:
            type: string
            enum: [auto, fraction, percent]
            description: Input viability scale (default auto).
          dose_unit:
            type: string
            description: Concentration unit label for plots and metadata (default uM).
      approval: write
      output: ./schemas/dose-response-result.json
      timeoutSeconds: 600
    - name: synergy
      description: >-
        Analyze 2D drug-combination dose-response matrices using HSA, Bliss
        Independence, Loewe Additivity, and ZIP (Zero Interaction Potency)
        synergy models (reported on the standard -100..+100 percent excess
        inhibition scale where >10 is synergistic and <-10 is antagonistic).
      run: [python, ./scripts/pgx_synergy.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Long-format 2D combination TSV/CSV with dose_a, dose_b, and viability columns (including dose_a=0 and dose_b=0 monotherapy/vehicle margins).
          output_dir:
            type: string
            format: project-path
            description: Directory for synergy summary, cell-level delta grid, and 2-panel synergy heatmap.
          dose_a_col:
            type: string
            description: Column for Drug A concentration (default dose_a).
          dose_b_col:
            type: string
            description: Column for Drug B concentration (default dose_b).
          viability_col:
            type: string
            description: Column for observed cell viability (default viability; auto-detects 0..1 vs 0..100%).
          drug_a_name:
            type: string
            description: Display name for Drug A (inferred from drug_a column if present, else DrugA).
          drug_b_name:
            type: string
            description: Display name for Drug B (inferred from drug_b column if present, else DrugB).
          sample_col:
            type: string
            description: Optional cell line / sample column when multiple combination matrices are stacked.
          viability_scale:
            type: string
            enum: [auto, fraction, percent]
            description: Input viability scale (default auto).
          dose_unit:
            type: string
            description: Concentration unit label (default uM).
      approval: write
      output: ./schemas/synergy-result.json
      timeoutSeconds: 600
    - name: biomarker
      description: >-
        Screen somatic alterations (MAF file or binary gene x sample matrix via
        Welch t-test, Mann-Whitney U, and Cohen's d) or continuous expression /
        pathway matrices (Spearman rho) against drug sensitivity metrics
        (IC50/LN_IC50/AUC/AAC) with explicit sensitivity polarity and BH FDR.
      run: [python, ./scripts/pgx_biomarker.py]
      args:
        type: object
        additionalProperties: false
        required: [drug_response, features, output_dir]
        properties:
          drug_response:
            type: string
            format: input-path
            description: Drug response matrix (samples x drugs or drugs x samples) or long table (sample, drug, metric).
          features:
            type: string
            format: input-path
            description: MAF somatic mutation table, binary alteration matrix (features x samples), or continuous expression/score matrix (features x samples).
          output_dir:
            type: string
            format: project-path
            description: Directory for biomarker association table and volcano/boxplot figure.
          feature_type:
            type: string
            enum: [auto, maf, binary, continuous]
            description: Genomic/molecular feature representation (default auto).
          metric_col:
            type: string
            description: Metric column in long-format drug_response (e.g. ln_ic50, ic50, auc, aac).
          sample_col:
            type: string
            description: Sample column in long-format drug_response (default sample).
          drug_col:
            type: string
            description: Drug column in long-format drug_response (default drug or compound).
          higher_is_sensitive:
            type: boolean
            description: Set true when higher drug metric values mean greater sensitivity (e.g. AAC or percent inhibition); leave false (default) for IC50, LN_IC50, or viability AUC.
          min_altered:
            type: integer
            minimum: 2
            description: Minimum altered and wild-type samples required per binary feature-drug test (default 3).
          fdr:
            type: number
            description: Benjamini-Hochberg FDR significance threshold (default 0.05).
      approval: write
      output: ./schemas/biomarker-result.json
      timeoutSeconds: 600
    - name: hts_qc
      description: >-
        Evaluate high-throughput drug screening plates (96/384-well) with
        Zhang's Z'-factor, SSMD, signal-to-background (S/B), control CV%,
        row/column ANOVA and edge-effect diagnostics, optional Tukey two-way
        median polish B-score normalization, percent inhibition, and robust MAD
        z_robust hit calling.
      run: [python, ./scripts/pgx_hts_qc.py]
      args:
        type: object
        additionalProperties: false
        required: [plates, output_dir]
        properties:
          plates:
            type: string
            format: input-path
            description: Plate screening TSV/CSV with plate (optional), well (e.g. A01..P24) or row+col, signal, and well_type (pos_ctrl, neg_ctrl, sample, blank).
          output_dir:
            type: string
            format: project-path
            description: Directory for hts_plate_qc_summary.tsv, hts_well_scored.tsv, and spatial QC heatmap.
          apply_bscore:
            type: boolean
            description: Apply Tukey two-way median polish B-score spatial bias correction for hit calling (default false).
          zscore_hit_cutoff:
            type: number
            description: Robust MAD |z_robust| cutoff for primary hit calling on sample wells (default 3.0).
      approval: write
      output: ./schemas/hts-qc-result.json
      timeoutSeconds: 600
    - name: predict_cohort
      description: >-
        Translate reference cell-line drug sensitivity models to clinical tumor
        expression cohorts using cross-platform gene-wise Z-score /
        empirical-Bayes harmonization and 5-fold cross-validated RidgeCV,
        flagging drugs with CV Pearson r < min_cv_r as low_predictability_warning.
      run: [python, ./scripts/pgx_predict_cohort.py]
      args:
        type: object
        additionalProperties: false
        required: [train_expr, train_response, clinical_expr, output_dir]
        properties:
          train_expr:
            type: string
            format: input-path
            description: Reference cell-line expression matrix TSV/CSV (genes x cell_lines or cell_lines x genes).
          train_response:
            type: string
            format: input-path
            description: Reference cell-line drug response matrix or long table (sample, drug, response).
          clinical_expr:
            type: string
            format: input-path
            description: Target clinical patient cohort expression matrix TSV/CSV.
          output_dir:
            type: string
            format: project-path
            description: Directory for cohort_predicted_response.tsv, drug_model_cv_metrics.tsv, and overview plot.
          min_cv_r:
            type: number
            description: Minimum 5-fold CV Pearson r on training cell lines to mark a drug model reliable (default 0.25).
          response_polarity:
            type: string
            enum: [lower_is_sensitive, higher_is_sensitive]
            description: Drug response polarity (lower_is_sensitive for IC50/AUC; higher_is_sensitive for AAC; default lower_is_sensitive).
      approval: write
      output: ./schemas/predict-cohort-result.json
      timeoutSeconds: 900
---

# Pharmacogenomics & Drug Screening Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for equations, censoring conventions, and drug-metric polarity rules.

## Workflow Chain

```
Raw HTS plate readout           ──pgx_hts_qc─────────▶ Z'-factor, SSMD, edge bias, B-score & primary hits
Raw plate / dose-viability data ──pgx_dose_response──▶ 4PL Hill fits (EC50, IC50, AUC, AAC, beyond_tested_range)
2D combination matrix           ──pgx_synergy────────▶ ZIP / Bliss / Loewe / HSA delta grids & peak synergy window
Drug IC50/AUC + MAF/Expression  ──pgx_biomarker──────▶ Sensitivity vs resistance biomarkers (Cohen's d / rho + FDR)
Cell-line train + tumor expr    ──pgx_predict_cohort─▶ Cross-platform RidgeCV patient sensitivity & CV reliability
```

## Hard Red Lines

1. **Do not report uncensored extrapolated $\text{IC}_{50}$ values**: if `beyond_tested_range` is `true`, report `> max_dose` and rank cell lines by `AUC` or `AAC`.
2. **Never invert sensitivity vs resistance polarity**: for `IC50`, `LN_IC50`, and viability `AUC`, negative $\Delta$ (altered < WT) means **collateral sensitivity** and positive $\Delta$ means **resistance**; for `AAC`, pass `--higher_is_sensitive` so polarity is standardized automatically.
3. **Require monotherapy margins in 2D synergy screens**: `pgx_synergy` requires `dose_a == 0` and `dose_b == 0` rows to anchor single-agent dose-response baselines.
