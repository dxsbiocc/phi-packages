---
name: grant-proposal
description: Draft, structure, and review competitive research grant proposals across China NSFC (青年C类, 面上, 重点; 机制探索型 & 方法算法构建型), NIH (R01, R21, K99/R00), NSF (Standard, CAREER), and ERC (Starting, Consolidator, Advanced). Enforces agency-specific rubrics, non-dependent Aims, risk-mitigation & alternative strategies, NSFC at-most-400-char abstract verification, and zero fabricated preliminary data or citations. Use for 国自然, NSFC标书, 基金申请书, 青年基金, 面上项目, NIH R01, R21, NSF grant, ERC grant, Specific Aims, and grant review. Not for journal manuscripts (use manuscript-writing) or detailed power/batch design (use study-design).
license: MIT
metadata: {"version": "1.1.0", "category": "文献与基金", "license": "MIT"}
---

# Research Grant Proposal Writing & Review (NSFC / NIH / NSF / ERC)

Design, draft, and review high-scoring grant proposals across **China NSFC** and international agencies (**NIH**, **NSF**, **ERC**) aligned with each agency's exact evaluation criteria, page/character limits, and panel psychology—**never fabricating preliminary data, sample sizes, or citations**.

## 1. Scope & Boundary

- **Use for**: Framing core scientific questions/hypotheses, structuring and drafting grant sections, and running mock panel reviews for **NSFC** (青年 C 类 / 面上 / 重点), **NIH** (R01 / R21 / K99-R00), **NSF** (Standard / CAREER), and **ERC** (StG / CoG / AdG).
- **Not for (Handoff)**:
  - Writing journal manuscripts -> `manuscript-writing`.
  - Deep hypothesis-competition or pre-data-collection power/batch design -> `study-design`.
  - Systematic literature synthesis -> `literature-review`.
  - Pure prose polishing or de-AI rewriting -> `manuscript-polishing` / `humanizer`.

## 2. Epistemic Contract & Stop Gates

1. **Zero Fabricated Preliminary Data or Track Record**: Never invent pilot figures, effect sizes, cohort `n`, prior publications, or equipment. Tag missing preliminary evidence `[PRELIMINARY DATA NEEDED / 待补充前期数据: <exact item>]` and unverified citations `[CITE / 待核验文献: <topic>]`.
2. **Agency & Mechanism Stop Gate**: If the user has not specified the **target agency & tier** (e.g., `NSFC 面上` vs `NSFC 青年 C 类` vs `NIH R01` vs `NSF CAREER` vs `ERC StG`) or core idea, **STOP and confirm** before drafting—agency rubrics and limits are mutually incompatible.
3. **Hard Character/Page Limit Gate**: For **NSFC Chinese abstracts (`<= 400` chars incl. punctuation)**, target **350–385 chars** and verify exact length via Python `len()` when available. Never exceed agency page budgets (`NIH 1-page Specific Aims`, `ERC 5-page B1`).

## 3. Phi Runtime & Context Routing

- **Workspace Grounding**: Inspect workspace CVs, biosketches, preliminary figures/tables, or `*.phi-artifact.json` sidecars to ground preliminary work and representative publications in verified artifacts.
- **On-Demand Agency References** (load ONLY the agency matching the user's proposal):
  - **China NSFC (青年 C 类 / 面上 / 重点 — 机制探索型 & 方法/算法构建型)** -> read [references/nsfc.md](references/nsfc.md).
  - **NIH (R01, R21, K99/R00)** -> read [references/nih-grants.md](references/nih-grants.md).
  - **NSF (Standard Research, CAREER)** -> read [references/nsf-grants.md](references/nsf-grants.md).
  - **ERC (Starting, Consolidator, Advanced)** -> read [references/erc-grants.md](references/erc-grants.md).

## 4. Core Workflow

1. **Agency, Tier & Track Lock**: Confirm agency, tier, track (if NSFC), and page/char limits; load ONLY the matching agency reference from Section 3.
2. **Logic Spine Blueprint**: Align `Problem & Gap -> Preliminary Rationale -> Central Hypothesis / Key Scientific Questions -> 2–4 Independent Synergistic Aims -> Pitfalls & Alternative Strategies -> Payoff`.
3. **Drafting / Mock Panel Review**: Draft section-by-section respecting agency rules, or score an existing proposal against the agency's official criteria, surfacing fatal panel red flags first.
4. **Verification**: Confirm zero fabricated preliminary data/citations, Aim independence, and exact abstract character/word limits.

## 5. Scientific Redlines

1. **No "Domino Aims"**: Aim 2/3 must remain executable and informative even if Aim 1's primary hypothesis turns out negative. Every Aim needs concrete **Pitfalls & Alternative Strategies**.
2. **Distinguish Scientific Questions from Technical Steps**: *"Build a knockout mouse"* or *"train a Transformer"* is a task, NOT a key scientific question.
3. **No Over-Claiming Completion**: Preliminary data must justify plausibility without implying the grant work is already finished.
4. **Zero Hype / AI Slogans**: Eliminate *"填补空白"*, *"开辟崭新道路"*, *"paradigm-shifting"*, and *"unprecedented"*.

## 6. Output Contract

1. **Grant Alignment Card** (Agency, Tier/Track, Review Attribute/Criteria, Core Hypothesis, Aim Independence Check).
2. **Drafted Section(s) or Scored Mock Panel Critique** (including `len()` check for NSFC abstracts).
3. **Gap & Compliance Manifest** (`[待补充前期数据]` / `[PRELIMINARY DATA NEEDED]`, `[CITE]` / `[待核验文献]`, agency compliance gaps).
