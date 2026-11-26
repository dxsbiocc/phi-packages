---
name: manuscript-writing
description: Draft or restructure scientific manuscripts across all journal tiers (Nature, Science, Cell, clinical journals, bioinformatics/engineering journals, and general IMRaD) from verified author-provided results, tables, and figures, plus initial-submission cover letters. Enforces hourglass narrative structure, 6-level causal claim-evidence calibration, reporting-guideline routing, and zero fabricated statistics or citations. Use for 论文写作, 章节起草, 论证重构, 初投 Cover Letter, abstract crafting, and IMRaD drafting. Not for line-editing existing prose (use manuscript-polishing), removing AI prose tells only (use humanizer), post-decision rebuttals (use peer-review-response), or statistical auditing (use statistical-review).
license: MIT
metadata: {"version": "1.1.0", "category": "论文写作", "license": "MIT"}
---

# Scientific Manuscript Writing (Universal High-Impact Standard)

Draft and restructure scientific manuscripts across any journal tier strictly grounded in **verified author-provided evidence**. Never invent experimental numbers, p-values, sample sizes, or citations.

## 1. Scope & Boundary

- **Use for**: Drafting new IMRaD sections (Title, Abstract, Introduction, Results, Discussion, Methods, Data/Code Availability) from verified tables/figures/reports, rebuilding manuscript argument flow, and writing initial-submission Cover Letters.
- **Not for (Handoff)**:
  - Line-level polishing, word-count trimming, or CN/EN translation of an existing draft -> `manuscript-polishing`.
  - Stripping AI prose tells only -> `humanizer`.
  - Pre-submission mock peer review -> `manuscript-review`.
  - Deep statistical auditing -> `statistical-review`.
  - Post-decision rebuttal & resubmission cover letter -> `peer-review-response`.
  - Converting a paper into slides -> `paper-to-slides`.

## 2. Epistemic Contract & Stop Gates

1. **Zero Fabrication**: Every number, effect size, CI, exact `n`, `P/FDR`, software version, and citation must trace to a user-provided artifact. Tag missing facts `[TBD: <exact item>]` and unverified references `[CITE: <topic>]`.
2. **Contradiction & Gap Stop Gate**: Build a **Claim-to-Evidence Map** before drafting prose. If source tables/figures contradict the user's draft claims OR >40% of core quantitative metrics are missing, **STOP at the Claim-to-Evidence Map** and list every conflict/`[TBD]` blocker before writing prose.

## 3. Phi Runtime & Context Routing

- **Artifact Context First**: Inspect workspace files and `.phi-artifact.json` (`scientificContext`: `valueSemantics`, `sampleIds`, `donorBySample`) so Results/Methods match exact upstream biological units.
- **Standalone**: Works without sibling skills or literature connectors.
- **On-Demand Branch References** (load ONLY the file needed for the user's request):
  - Targeting specific journal formats (*Nature/Science/Cell*, clinical, bioinformatics) or drafting an Initial Cover Letter -> read [references/journal-tiers.md](references/journal-tiers.md).
  - Drafting Methods / Compliance for RCT, observational, animal, ML/AI, diagnostic, systematic review, or MR studies -> read [references/reporting-guidelines.md](references/reporting-guidelines.md).

## 4. Core Workflow

1. **Intake**: Confirm target journal tier, study design, target section(s), and language (EN/ZH). Load only the relevant branch reference if applicable.
2. **Claim-to-Evidence Map**: Map `Subsection -> Core Claim -> Causal Level -> Supporting Artifact -> Exact Metrics -> Gaps/Conflicts`. Trigger Stop Gate if blocked.
3. **Drafting (Hourglass Architecture)**:
   - **Title & Abstract**: Broad problem -> specific gap -> concrete advance -> quantitative finding -> significance.
   - **Introduction (Funnel)**: Known landscape -> unresolved gap -> why it matters -> *"Here, we show..."* preview.
   - **Results**: One testable claim per subsection; lead with finding, support with effect size + `95% CI` (or `SD/IQR`) + biological `n` + `P/FDR`, and cite exact figure panels.
   - **Discussion (Reverse Funnel)**: Lead with conceptual advance (do not repeat Results) -> reconcile literature -> scope/implications -> candid limitations -> falsifiable next steps.
   - **Methods**: Reproducible protocol, exact software versions, statistical tests, multiple-testing corrections, and reporting-guideline items.
4. **Pre-Delivery Audit**: Verify causal-level verbs, biological vs technical `n`, and compile a consolidated `[TBD]` / `[CITE]` manifest.

## 5. Scientific Redlines

1. **6-Level Causal Ladder**:
   - L1 Descriptive (*observed, detected*), L2 Cross-sectional (*associates with, correlates with, is enriched in*), L3 Longitudinal/Predictive (*predicts, stratifies*), L4 Quasi-experimental/MR (*is consistent with a causal effect of*), L5 Interventional/RCT (*increases, reduces*), L6 Perturbation + Rescue (*is required for, mediates*).
   - **Hidden Causal Verb Blacklist** (forbidden for L1–L3): *drives, promotes, triggers, mediates, protects against, reverses, impairs, orchestrates*.
2. **Pseudoreplication Guard**: Never report cell/ROI counts as biological sample size; report both biological units (`n = 12 donors`) and technical observations.
3. **No AI Prose Tells**: Avoid *"In recent years..."*, *"plays a pivotal role"*, forced triads, and theatrical closers.

## 6. Output Contract

1. **Claim-to-Evidence Map** (table).
2. **Drafted Section(s)** (Markdown or LaTeX matching project format, exact panel citations).
3. **Compliance & Epistemic Manifest** (guideline applied + consolidated `[TBD]` and `[CITE]` list).
