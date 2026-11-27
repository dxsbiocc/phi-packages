---
name: research-provenance
description: "Manage multi-step scientific investigation state, falsifiable hypothesis tracking, end-to-end data-to-claim lineage, cross-literature and cross-modality contradiction matrices, and formal stopping criteria across Phi workflows. Maintains lightweight research_ledger.json and PROVENANCE.md artifacts mapping input SHA-256 hashes, tools, parameters, artifacts, and effect sizes to manuscript claims to prevent orphan numbers or silent goal-post shifting. Use for 研究溯源, 假说台账, 证据链审计, 矛盾调和, research provenance, hypothesis ledger, evidence lineage, and contradiction matrix. Not for initial study design (use study-design) or standalone paper critique (use literature-review)."
license: MIT
metadata: {"version": "1.0.0", "category": "科学分析", "license": "MIT"}
---

# Scientific Investigation Provenance & Hypothesis Ledger

Maintain auditable, reproducible epistemic state across multi-step, multi-omics, and literature-integrated research campaigns—preventing orphan numbers, post-hoc goal-post shifting, and unresolved cross-study or cross-modality contradictions.

## 1. Scope & Boundary

- **Use for**:
  - **Hypothesis & Belief-Update Tracking**: Maintaining `H0` + competing `H1..Hk` hypotheses with pre-registered falsification thresholds, posterior statuses, and explicit evidence-update logs.
  - **End-to-End Data-to-Claim Lineage**: Linking every quantitative claim to its `input_sha256`, Phi tool/wrapper, exact parameters, generated artifact path, and `(effect_size, CI, P/FDR)`.
  - **Cross-Literature & Cross-Modality Contradiction Auditing**: Systematically adjudicating disagreements (e.g., RNA vs protein, *in vitro* vs patient cohort, Paper A vs Paper B) by power, batch/confounding, assay modality, and bias risk.
  - **Investigation Stopping Criteria**: Evaluating whether evidence has converged or requires targeted follow-up before moving to manuscript drafting.
- **Not for (Handoff)**:
  - Pre-data-collection experimental design or power calculation -> `study-design`.
  - Deep single-paper critical reading or systematic literature search -> `literature-review`.
  - Statistical assumption auditing of a single manuscript -> `statistical-review`.
  - Drafting final IMRaD prose -> `manuscript-writing`.

## 2. Epistemic Invariants (Hard Gates)

1. **Zero Orphan Numbers**: Any statistic (`FC`, `HR`, `OR`, `AUC`, `Z'`, `IC50`, `P/FDR`) promoted to a conclusion or manuscript claim MUST have a complete `lineage` entry (`input_sha256 -> tool_or_wrapper -> exact_parameters -> artifact_path -> effect_size_ci_p -> manuscript_claim_id`). Unlinked numbers are tagged `[ORPHAN CLAIM — UNVERIFIED]`.
2. **Pre-Registered Falsification & No Silent Edits**: Once a hypothesis (`H1..Hk`) is logged with a falsification threshold, never silently rewrite the threshold after seeing negative or messy results. Log an explicit `evidence_updates` entry and transition `status` honestly (`supported`, `refuted`, `inconclusive`, or `confounded`).
3. **Hash-Verified Freshness**: When upstream raw inputs or scripts change (`sha256` mismatch against `research_ledger.json`), flag all downstream claims derived from the old hash as `[STALE LINEAGE — RE-RUN REQUIRED]`.
4. **Explicit Contradiction Adjudication**: Never cherry-pick supportive evidence while ignoring discordant modalities or conflicting papers. Every known conflict must enter the Contradiction Matrix with an explicit adjudication verdict.

## 3. On-Demand Reference

- For the exact `research_ledger.json` schema template, `PROVENANCE.md` layout, and a worked multi-plugin example (`literature-library` + `single-cell` + `wetlab-assay` + `pharmacogenomics`) -> read [references/ledger-and-contradiction-schema.md](references/ledger-and-contradiction-schema.md).

## 4. Core Workflow

1. **Initialize or Load Ledger**: Create or inspect `research_ledger.json` and `PROVENANCE.md` in the investigation directory. Register `H0` (batch/artifact/reverse-causality null) plus `H1..Hk` with pre-registered falsification criteria and explicit stopping rules.
2. **Append Verified Lineage Records**: After each Phi tool/wrapper run (e.g., `sc_de_markers`, `bulk_deseq2`, `pgx_hts_qc`, `wa_qpcr_quant`, `wa_immunoassay_curve`, `lit_get`), record `input_sha256`, exact CLI/tool arguments, output artifact path (plus `.phi-artifact.json` sidecar if emitted), and observed `effect_size_ci_p`.
3. **Update Hypotheses & Reconcile Contradictions**: Update each hypothesis's posterior state with quantitative/qualitative evidence weights (`log_likelihood_ratio` or structured strength tier) and populate the Contradiction Matrix whenever two datasets, assays, or publications disagree.
4. **Evaluate Stopping Criteria Gate**: Before declaring an investigation complete or handing off to `manuscript-writing`, check the formal stopping rules:
   - **CONVERGE**: Primary `H` crosses pre-registered confirmation threshold across `>= 2` orthogonal lines/modalities, `H0` confounders are ruled out, all high-severity contradictions are adjudicated, and `0` orphan claims remain.
   - **CONTINUE / BRANCH**: Evidence is `inconclusive` or `confounded` with an actionable discriminating test remaining.
   - **STOP (FALSIFIED / UNDERPOWERED)**: Primary `H` meets falsification criteria or remaining experiments exceed feasible power/budget—report negative/inconclusive result cleanly.

## 5. Output Contract

Every `research-provenance` update or audit produces:
1. **`research_ledger.json`** (machine-readable ledger following the schema in `references/ledger-and-contradiction-schema.md`).
2. **`PROVENANCE.md`** (human-readable dashboard containing: Hypothesis Status Board, Data-to-Claim Lineage Table, Contradiction Matrix, and Stopping-Gate Verdict).
3. **Audit Alert List** (explicitly flagging any `[ORPHAN CLAIM]`, `[STALE LINEAGE]`, or `[UNRESOLVED CONTRADICTION]` blockers).
