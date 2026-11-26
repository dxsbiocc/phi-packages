---
name: literature-review
description: Conduct systematic or scoping literature synthesis across >= 5 studies, build structured evidence matrices, reconcile conflicting findings, verify retraction/erratum and preprint-to-published status, and map actionable research gaps. Routes retrieval through Phi literature tools (lit_search, lit_query, lit_curate) or database MCPs and enforces a strict No-Retrieval Stop Gate to prevent fabricated citations. Use for 文献综述, 文献调研, 证据矩阵, systematic review, scoping review, and research gap analysis. Not for deep reading of 1-3 papers (use paper-reading) or full IMRaD drafting (use manuscript-writing).
license: MIT
metadata: {"version": "1.1.0", "category": "文献与基金", "license": "MIT"}
---

# Systematic & Scoping Literature Review (>= 5 Studies)

Synthesize multi-study evidence (`>= 5` papers) into a rigorous, verifiable **Evidence Matrix** and thematic narrative—never a chronological laundry list and **never a single fabricated citation**.

## 1. Scope & Boundary

- **Use for**: Scoping reviews, thematic syntheses, and PRISMA-guided systematic reviews across `>= 5` studies; comparative Evidence Matrices; contradiction reconciliation; research gap analysis.
- **Not for (Handoff)**:
  - Deep figure-by-figure reading of `1–3` papers -> `paper-reading`.
  - Drafting a full original research manuscript -> `manuscript-writing`.
  - Writing an NSFC grant proposal -> `nsfc-proposal`.

## 2. Epistemic Contract & Stop Gates

1. **No-Retrieval Stop Gate (Hard Stop)**:
   - **Never fabricate DOIs, PMIDs, authors, years, or study findings from parametric memory.**
   - If NO active retrieval tool (`lit_search` / `pubmed` / `biorxiv` / `europepmc`) is available AND the user provides NO local PDFs/bibliography, **STOP immediately**. Output only a **Search Strategy Blueprint** (PICO/PECO Boolean queries, inclusion/exclusion criteria, screening schema) and instruct the user to enable a literature connector or attach PDFs.
2. **Primary-Source Citation Discipline**: Never cite a secondary review article for a primary experimental number or cohort size. Trace quantitative claims to the primary study, or tag `[SECONDARY REVIEW — VERIFY PRIMARY SOURCE]`.
3. **Retraction & Preprint Version Audit**: Check retraction/erratum/expression-of-concern status and preprint-to-published resolution whenever tools allow. Quarantine retracted studies and tag unresolved preprints `[PREPRINT]`.

## 3. Phi Runtime & Context Routing

Follow this strict retrieval priority chain:
1. **`literature-library` plugin (if enabled)**: `lit_search` / `lit_query` for retrieval, `lit_get(with_fulltext=True)` for key papers, `lit_curate(action="versions")` for preprint -> published resolution.
2. **Database MCPs (if enabled)**: `pubmed`, `biorxiv`, `europepmc`, `openalex` (log exact query string, date, and hit counts).
3. **Local Workspace Corpus**: User-supplied PDFs, `.bib`, or CSV/Zotero exports (`Scoped to user-supplied corpus, k = ...`).
4. **Stop Gate**: If 1–3 are absent -> halt and deliver the Search Strategy Blueprint only.

## 4. Core Workflow & Evidence Matrix Schema

1. **Protocol**: Frame scope via **PICO/PECO** or **Method-Benchmark-Domain**; define inclusion/exclusion criteria and date range.
2. **Retrieval & Integrity Screening**: Run reproducible queries; check retractions/errata; deduplicate preprints vs published versions; separate primary studies from reviews.
3. **Build Evidence Matrix**:
   `| Study (Author, Year, ID) | Design & Causal Level | System/Cohort (n) | Primary Metric (Effect Size, CI, P) | Key Limitation | Status |`
4. **Thematic Synthesis**: Group by mechanism, methodology, or clinical axis (never *"Author A found X. Author B found Y."*). Explicitly explain **contradictions/heterogeneity** (stage, assay, sample size, statistical power) and **publication bias**.
5. **Falsifiable Research Gaps**: Conclude with 2–4 specific gaps directly grounded in the matrix.

## 5. Worked Micro-Example (Thematic Reconciliation vs Laundry List)

- *Weak (Laundry list + secondary review)*: *"Smith (2021) found SPP1+ macrophages predicted poor survival. Wang (2023) found no effect. Review Brown (2024) says they cause resistance in 70% of tumors."*
- *Strong (Thematic synthesis)*: *"Across four GI tumor cohorts (`k = 4`, `n = 518`), elevated *SPP1*+ macrophage abundance associated with shorter OS in three studies (`HR` `1.62–2.14`; Smith 2021 `PMID:...`, Jones 2022 `PMID:...`, Park 2023 `PMID:...`), whereas Wang 2023 (`PMID:...`, `n = 46`) reported a null result (`HR = 1.12`, `95% CI [0.71, 1.77]`), likely reflecting Stage I–II enrollment and low event power (`19` deaths). Secondary review Brown 2024 was excluded from quantitative pooling."*

## 6. Output Contract

1. **Search & Screening Log** (sources, query/date, `k` included, retracted/preprint-resolved notes).
2. **Structured Evidence Matrix** (every study carries a verified `DOI`/`PMID`/preprint ID or local filename).
3. **Thematic Synthesis & Contradiction Analysis** + **2–4 Falsifiable Research Gaps**.
