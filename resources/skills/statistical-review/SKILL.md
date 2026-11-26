---
name: statistical-review
description: Audit statistical methods, study design, experimental units, survival/omics/ML assumptions, uncertainty reporting, and figure statistics in scientific manuscripts, analysis notebooks, or reviewer responses. Catches pseudoreplication, PH/competing-risk/immortal-time bias, compositional and scRNA-seq double-dipping errors, data leakage, and p-value misuse. Use for 统计审查, 统计方法核查, 图注统计, biostatistics audit, and reviewer statistical defense. Not for general manuscript prose editing (use manuscript-polishing) or full editorial novelty review (use manuscript-review).
license: MIT
metadata: {"version": "1.1.0", "category": "科学分析", "license": "MIT"}
---

# Scientific Statistical & Methodological Audit

Audit the statistical validity, assumption compliance, and reporting rigor of scientific manuscripts, figure legends, and analysis scripts. **Never fabricate calculations or eyeball p-values from plots.**

## 1. Scope & Boundary

- **Use for**: Auditing experimental units (`biological` vs `technical`), test selection, assumption checks, multiple-testing corrections, uncertainty reporting, domain-specific statistical redlines, and writing rigorous Statistical Methods / figure legends.
- **Not for (Handoff)**:
  - Pre-data-collection hypothesis framing, batch design, or power planning -> `study-design`.
  - Broad editorial/novelty mock peer review -> `manuscript-review`.
  - Drafting IMRaD narratives -> `manuscript-writing`.
  - Polishing prose -> `manuscript-polishing`.
  - Full rebuttal letter formatting -> `peer-review-response`.

## 2. Epistemic Contract & Stop Gates

1. **No Eyeballing Plots**: Never visually estimate or "recompute" p-values, HRs, or error bars from bar charts or Kaplan-Meier curves.
2. **Unverifiable-from-Text Gate**: If a manuscript omits test assumptions, independence details, missing-data handling, or code, label the item `[UNVERIFIABLE FROM TEXT — INSPECT CODE/RAW DATA]`—never assume it passes or fails without evidence.
3. **Compute Only When Data + Request Exist**: Audit text/code by default; run new calculations only when the user supplies raw data/scripts and explicitly requests re-analysis.

## 3. Phi Runtime & Context Routing

- **Artifact & Code Grounding**: Inspect upstream tool outputs (`scientificContext`: `valueSemantics`, `sampleIds`, `donorBySample`, `batch`), `*.phi-artifact.json` sidecars (`provenance`, `parameters`), and workspace scripts (`*.R`, `*.py`, `*.ipynb`) when available.
- **On-Demand Domain References** (load ONLY the branch matching the user's study):
  - **Survival / Clinical / Epidemiology** (Cox, KM, competing risks, immortal time, EPV, subgroups, confounding) -> read [references/survival-clinical.md](references/survival-clinical.md).
  - **Bulk RNA-seq / scRNA-seq / Spatial / Microbiome / Proteomics** (pseudoreplication, raw counts vs TPM, compositional data, double dipping, LOD/MNAR) -> read [references/omics-singlecell.md](references/omics-singlecell.md).
  - **Machine Learning / Clinical Prediction / AI** (data leakage, nested CV, calibration, class imbalance) -> read [references/ml-prediction.md](references/ml-prediction.md).

## 4. Core Workflow

1. **Experimental Unit Map**: Separate biological units (`patient`, `donor`, `mouse`, `cage`) from technical observations (`cells`, `ROIs`, `wells`). Flag any mismatch immediately.
2. **Universal + Domain Redline Audit**:
   - Check universal rules: exact `n`, effect size + `95% CI` (or `SD/IQR`), exact `P/FDR`, two-sidedness, multiple-testing control (`BH-FDR` / `Bonferroni`), no bar-chart-only continuous data (`n <= 25`), no `"significant in A, NS in B"` without an interaction test.
   - Load ONLY the relevant domain reference from Section 3 and run its domain-specific checks.
3. **Prioritized Audit Report**: Classify each finding as `Critical (Invalidates Inference)`, `Major (Needs Re-analysis / Sensitivity Check)`, or `Minor (Reporting Clarification)` with a concrete replacement method.

## 5. Output Contract

1. **Experimental Unit & Data-Semantics Summary** (biological `n` vs technical observations).
2. **Statistical Audit Table** (`ID | Severity | Location | Flaw | Verification Status | Recommended Fix`).
3. **Drop-In Revised Statistical Methods / Figure Legend** (with `[TBD]` for unverified numbers).
4. **Unverifiable Items List** (`[UNVERIFIABLE FROM TEXT — INSPECT CODE/RAW DATA]`).
