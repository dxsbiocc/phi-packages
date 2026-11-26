---
name: study-design
description: Formulate falsifiable competing scientific hypotheses (strong inference) and design bias-resistant pre-data-collection experiments and cohorts across wet-lab, preclinical animal, clinical, and multi-omics settings. Covers crucial discriminating experiments, randomization, blinding, blocking, batch-treatment de-confounding, control architecture, and sample-size/power sensitivity planning. Use for 实验设计, 科学假说, 选题构思, 样本量估算, 效能分析, 随机化与区组设计, hypothesis generation, experimental design, and power analysis. Not for post-hoc statistical auditing of existing results (use statistical-review) or full grant drafting (use grant-proposal).
license: MIT
metadata: {"version": "1.1.0", "category": "科学分析", "license": "MIT"}
---

# Scientific Hypothesis Generation & Experimental Design

Bridge the pre-data-collection phase end-to-end: turn preliminary observations into **falsifiable competing hypotheses** and **bias-resistant, properly powered experimental designs** before a single sample is collected.

## 1. Scope & Boundary

- **Use for**:
  - **Phase A (Hypothesis & Strong Inference)**: Framing competing hypotheses (`H0` artifact/null vs `H1`/`H2` mechanisms), causal DAGs, novelty gating, and crucial discriminating experiments.
  - **Phase B (Experimental & Power Design)**: Biological vs technical experimental units, control architecture, randomization/blinding/blocking, batch-treatment de-confounding, and sample-size/power sensitivity planning.
- **Not for (Handoff)**:
  - Post-hoc statistical auditing of already-collected data or manuscripts -> `statistical-review`.
  - Literature reading or systematic synthesis -> `literature-review`.
  - Full grant drafting (NSFC / NIH / NSF / ERC) -> `grant-proposal`.
  - Manuscript Methods/Results drafting -> `manuscript-writing`.

## 2. Epistemic Contract & Stop Gates

1. **Zero Pseudo-Hypotheses**: Never output unfalsifiable vagueness (*"Gene X plays a complex role"*). Every hypothesis must state **Mechanism -> Direction -> Observable Prediction -> Falsification Criterion**, and screen at least 3 null/artifact explanations (`H0`: batch, cell-composition shift, reverse causality/confounding).
2. **No Fabricated Power Numbers**: Never invent a magic sample size (`"use n = 6/group"`). Every power estimate must state **primary endpoint**, **minimum biologically meaningful effect size**, **expected variance / event rate**, **two-sided alpha**, and **target power (`>= 80%`)**. If variance/effect priors are missing, output a **sensitivity grid** and tag `[PILOT ESTIMATE REQUIRED]`.
3. **Confounding Hard Stop**: If a proposed schedule runs all Treatment samples in Batch 1 and all Control samples in Batch 2, **halt and flag an unrecoverable confounding error**—no downstream batch correction can separate biology from batch.

## 3. Phi Runtime & Context Routing

- **Workspace Grounding**: Inspect workspace preliminary DE tables, pilot outputs, or `*.phi-artifact.json` sidecars when available; if Python/R is active and explicit priors are supplied, offer to compute exact power curves.
- **On-Demand Design References** (load ONLY the branch matching the user's experimental system):
  - **Wet-Lab, Cell & Preclinical Animal Studies** (ARRIVE 2.0, cage/litter blocking, factorial/rescue controls) -> read [references/wetlab-preclinical.md](references/wetlab-preclinical.md).
  - **Clinical Cohorts & High-Throughput Omics** (sequencing batch balance, donor-level power, stratification, survival events) -> read [references/clinical-omics.md](references/clinical-omics.md).

## 4. Core Workflow

1. **Hypothesis & Strong Inference (`H0` vs `H1` vs `H2`)**: Distill the core anomaly, sketch a minimal **causal DAG** (confounders, mediators, colliders), screen artifact/null explanations (`H0`), formulate 2–3 competing mechanisms, and define the **crucial discriminating experiment** whose outcomes separate them cleanly.
2. **Estimand & Experimental Unit Lock**: Define primary endpoint (continuous / binary / time-to-event) and true **biological unit** (`animal`, `cage`, `donor`, `patient`, `independent culture`)—never confuse cells, wells, or histological fields with `n`.
3. **Controls, Randomization, Blinding & Blocking**: Specify vehicle/negative, positive, and rescue controls; balance arms across sex (SAGER/SABV), litter/age, operator, plate, and sequencing run/lane.
4. **Power & Pre-Specification**: Build a sample-size sensitivity table (with 10%–20% attrition buffer), pre-specify primary statistical model, and lock QC/exclusion rules *before* data collection.

## 5. Scientific Redlines

1. **Biological Unit = Randomization Unit**: Cage-level dosing/co-housing makes the **cage** the experimental unit, not individual mice.
2. **Never Confound Batch with Group**: Every plate, dissection day, or sequencing lane must contain balanced arms.
3. **Pre-Specify Primary Endpoint & Exclusions**: Lock primary outcome and QC rules *a priori* to prevent p-hacking.
4. **Sex-as-a-Biological-Variable (SABV / SAGER)**: Include both sexes as a blocking factor unless scientifically unjustified.

## 6. Output Contract

1. **Competing Hypotheses & Crucial Experiment Table** (when framing hypotheses: `H0/H1/H2 | Prediction | Falsification Criterion`).
2. **Study Design Summary Card** (Primary Endpoint, Biological Unit, Controls, Randomization/Blinding/Blocking, Batch Layout).
3. **Sample-Size & Power Sensitivity Table** (explicit assumptions, attrition buffer, pre-specified primary test).
4. **Fatal-Flaw Prevention Checklist** (key confounding/pseudoreplication risks mitigated).
