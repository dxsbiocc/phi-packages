---
name: manuscript-review
description: Conduct comprehensive pre-submission manuscript readiness audits including mock peer review (Handling Editor + Domain Reviewers), claim-to-citation faithfulness & BibTeX verification, reporting-guideline & image/reagent integrity checks (CONSORT 2025, STROBE, ARRIVE 2.0, TRIPOD+AI, SAGER), and FAIR data/code availability compliance. Use for 模拟审稿, 投稿前预审, 引用核查, BibTeX核验, 数据可用性声明, 代码开源声明, manuscript critique, citation audit, and FAIR compliance. Not for deep statistical auditing (use statistical-review), answering received reviewer comments (use peer-review-response), or rewriting prose (use manuscript-polishing).
license: MIT
metadata: {"version": "1.1.0", "category": "论文写作", "license": "MIT"}
---

# Pre-Submission Manuscript Review, Citation Audit & FAIR Compliance

Run a complete pre-submission readiness gate—mock refereeing, reporting/reagent/image integrity, citation & BibTeX verification, and FAIR data/code availability compliance—*before* submitting to a journal.

## 1. Scope & Boundary

- **Use for**:
  - **Mock Peer Review** (Editor + Domain Reviewers): Novelty, controls, 6-level causal calibration, reporting-guideline & reagent/image integrity, and desk-reject risk triage.
  - **Citation & BibTeX Audit**: Claim-citation faithfulness, phantom/hallucinated DOIs/PMIDs, retractions/errata, preprint-to-published upgrades, and `.tex`/`.bib` sync.
  - **Data & Code Availability (FAIR)**: Public vs controlled repository routing (`GEO`, `EGA`, `dbGaP`, `GSA-Human`, `PRIDE`, `Zenodo`), human-data privacy gates, and compliant availability statements.
- **Not for (Handoff)**:
  - Deep statistical/methodological auditing -> `statistical-review` (flag obvious red flags inline if standalone).
  - Answering received journal reviewer reports -> `peer-review-response`.
  - Drafting or polishing manuscript prose -> `manuscript-writing` / `manuscript-polishing`.
  - Reading/synthesizing external published papers -> `literature-review`.

## 2. Epistemic Contract & Stop Gates

1. **Evidence-Grounded Critique Only**: Anchor every critique to an exact section, panel, table, or quote. Never visually eyeball p-values from plots; mark missing assumptions `[UNVERIFIABLE FROM TEXT — INSPECT CODE/RAW DATA]`.
2. **Partial-Draft Degradation Gate**: If only an Abstract or outline is supplied, **STOP full referee scoring** and degrade to **Desk-Triage Screening** (`[FULL MANUSCRIPT REQUIRED]`).
3. **Zero Memory-Based Citation Verification**: Never mark a DOI/PMID/claim `"Verified"` from memory alone—external verification requires active `lit_*`/`pubmed`/`biorxiv` tools or local PDFs; otherwise run Structural Audit Mode only (`[UNVERIFIED — RETRIEVAL TOOL REQUIRED]`).
4. **Zero Fake Accessions & Human-Privacy Hard Stop**: Never invent real-looking accessions (`GSE######` must be `[ACCESSION TBD]`). Forbid open public deposition of identifiable human raw genomic reads (`FASTQ/BAM/VCF`) or PHI—route strictly to controlled repositories (`EGA`/`dbGaP`/`GSA-Human`).

## 3. Phi Runtime & Context Routing

- **Workspace Audit**: Inspect manuscript files, `.bib`, figures, environment locks, and `*.phi-artifact.json` sidecars (if present).
- **On-Demand Compliance References** (load ONLY if the user requests or the audit triggers these modules):
  - **Citation Faithfulness, Retraction & BibTeX Audit** -> read [references/citation-audit.md](references/citation-audit.md).
  - **FAIR Repository Routing, Human Privacy & Data/Code Availability Statements** -> read [references/data-code-availability.md](references/data-code-availability.md).

## 4. Core Workflow (Mock Peer Review)

1. **Desk-Triage & Guideline Match**: Confirm target journal tier and apply the study's standard (`CONSORT 2025` / `STROBE` / `ARRIVE 2.0` / `TRIPOD+AI` / `STARD` / `PRISMA` / `STROBE-MR`) plus universal integrity checks (`SAGER` sex/gender, STR/mycoplasma, RRIDs, uncropped WB full blots, scale bars, AI disclosure).
2. **Claim-to-Evidence & Causal Audit**: Check headline claims against the 6-level causal ladder (`L1` Descriptive, `L2` Cross-sectional, `L3` Longitudinal/Predictive, `L4` Quasi-experimental/MR, `L5` Interventional/RCT, `L6` Perturbation + Rescue); flag blacklisted causal verbs (*drives, promotes, triggers, mediates, reverses, impairs*) backed only by observational (`L1–L3`) data.
3. **Optional Citation / FAIR Modules**: Load `references/citation-audit.md` and/or `references/data-code-availability.md` when auditing references or availability statements.
4. **Prioritized Referee Report**: Put fatal flaws first. Every Major Comment (`M1`, `M2`, ...) must state: `Location | Rejection Risk | Critique | Minimum Viable Fix`.

## 5. Output Contract

1. **Handling Editor Summary** (journal fit, novelty, verdict: `Ready / Minor Revision / Major Revision / High Desk-Reject Risk`).
2. **Reporting & Integrity Checklist** (guideline, SAGER, STR/RRID, WB/images, Data/Code).
3. **Major Comments (`M1`, `M2`, ...)** & **Minor Comments (`m1`, `m2`, ...)** + `[UNVERIFIABLE FROM TEXT]` items.
4. **Citation Audit Table / Availability Statements** (when those modules were invoked).
