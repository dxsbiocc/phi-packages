---
name: literature-review
description: Perform both deep critical paper reading (1 to 4 papers, figure-by-figure evidence chain and appraisal) and multi-study systematic or scoping literature synthesis (5+ studies, evidence matrix, contradiction reconciliation, and gap analysis). Routes retrieval via Phi literature tools (lit_get, lit_search, lit_query, lit_curate) or PubMed/bioRxiv MCPs and enforces strict Abstract-Only and No-Retrieval Stop Gates plus retraction/preprint checks. Use for 论文精读, 读论文, 文献综述, 文献调研, 证据矩阵, paper critique, systematic review, and scoping review. Not for full original manuscript drafting (use manuscript-writing) or pre-submission refereeing of your own manuscript (use manuscript-review).
license: MIT
metadata: {"version": "1.1.0", "category": "文献与基金", "license": "MIT"}
---

# Critical Paper Reading & Systematic Literature Review

Cover the entire literature-appraisal spectrum under one epistemic roof: **deep critical reading of 1–3 papers** and **cross-study thematic synthesis of `>= 5` papers**—with zero fabricated citations or hallucinated figures.

## 1. Scope & Boundary

- **Use for**:
  - **Deep Reading Mode (`1–4` papers)**: Figure-by-figure evidence chain, causal calibration, methodological critique, overclaim detection, reproducibility audit, and head-to-head comparison.
  - **Synthesis Mode (`5+` papers)**: PICO/PECO scoping or PRISMA-guided systematic reviews, comparative Evidence Matrices, contradiction/heterogeneity reconciliation, and falsifiable research gaps.
- **Not for (Handoff)**:
  - Pre-submission mock peer review or BibTeX audit of the user's own manuscript -> `manuscript-review`.
  - Drafting an original IMRaD research paper -> `manuscript-writing`.
  - Framing new competing hypotheses & experimental plans -> `study-design`.
  - Converting a verified paper into slides -> `paper-to-slides`.

## 2. Epistemic Contract & Stop Gates

1. **No-Retrieval Stop Gate (Synthesis Mode)**: Never fabricate DOIs, PMIDs, authors, years, or study findings from parametric memory. If NO retrieval tool (`lit_search` / `pubmed` / `biorxiv`) is active AND NO local PDFs/bibliography are supplied, **STOP** and output only a **Search Strategy Blueprint** (PICO Boolean queries + inclusion/exclusion criteria).
2. **Abstract-Only Degradation Gate (Deep Reading Mode)**: Never fabricate Fig. 1–5 maps, sample sizes, or Methods critiques from an abstract alone. Without full text, **STOP full Reading Card generation** and output strictly an **Abstract-Only Triage Card** marked `[FULL TEXT REQUIRED]`.
3. **Primary-Source & Retraction Discipline**: Never cite a secondary review for a primary experimental number. Check retraction/erratum status and preprint-to-published resolution whenever tools allow; quarantine retracted articles and tag unresolved preprints `[PREPRINT — NOT PEER-REVIEWED]`.

## 3. Phi Runtime & Context Routing

Follow this strict retrieval priority chain:
1. **`literature-library` plugin (`LiteratureLibrary` agent, if enabled)**: `lit_search` (sources: `pubmed`, `europepmc`, `openalex`, `biorxiv`) / `lit_query` for discovery, `lit_get(..., with_fulltext=True)` for locally saved PDF full text, `lit_curate(action="versions", http=True)` for preprint -> published resolution.
2. **Database MCPs (if enabled)**: `pubmed`, `biorxiv` (or user-configured external literature MCPs; log exact query string and search date).
3. **Local Workspace Corpus**: User-supplied PDFs, Markdown, `.bib`, or CSV exports.
4. **Stop Gate Fallback**: Trigger the applicable Stop Gate from Section 2 if 1–3 fail.

**On-Demand Mode References** (load ONLY the branch matching paper count):
- **1–4 papers (Deep Critical Reading / Comparison)** -> read [references/single-paper-card.md](references/single-paper-card.md).
- **`5+` papers (Systematic / Scoping Synthesis & Evidence Matrix)** -> read [references/multi-paper-synthesis.md](references/multi-paper-synthesis.md).

## 4. Core Workflow

1. **Route by Scope & Verify Sources**: Determine Deep Reading (`1–4`) vs Synthesis (`5+`); verify full-text/retrieval availability and load ONLY the matching branch reference.
2. **Integrity & Version Screening**: Verify publication status (peer-reviewed vs preprint), check retractions/errata, and separate primary empirical studies from secondary reviews.
3. **Structured Appraisal / Synthesis**: Execute the exact template in the loaded branch reference—separating author claims from demonstrated evidence on the 6-level causal ladder (`L1` Descriptive, `L2` Cross-sectional, `L3` Longitudinal/Predictive, `L4` Quasi-experimental/MR, `L5` Interventional/RCT, `L6` Perturbation + Rescue).
4. **Actionable Transfer / Research Gaps**: Conclude with concrete reusable methods (Deep Reading) or 2–4 falsifiable knowledge gaps anchored to the Evidence Matrix (Synthesis).

## 5. Output Contract

- **Deep Reading (`1–3` papers)**: Citation & Status Header, Executive Summary, Figure-by-Figure Evidence Chain (or `[FULL TEXT REQUIRED]` Triage Card), Critical Appraisal (overclaims/controls/reproducibility), Actionable Takeaways.
- **Synthesis (`>= 5` papers)**: Search & Screening Log, Structured Evidence Matrix (verified IDs only), Thematic Synthesis & Contradiction Analysis, 2–4 Falsifiable Research Gaps (or Search Strategy Blueprint if stopped).
