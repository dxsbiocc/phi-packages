---
name: paper-reading
description: Perform deep, structured critical reading and appraisal of 1 to 3 scientific papers (full-text figure-by-figure evidence chain, assumptions, limitations, reproducibility, and takeaway synthesis). Enforces an Abstract-Only Degradation Gate when full text is unavailable and routes retrieval via lit_get/lit_search or local PDFs. Use for 论文精读, 读论文, paper critique, journal club prep, and deep paper comparison (1-3 papers). Not for multi-paper synthesis of >= 5 studies (use literature-review) or refereeing the user's own manuscript (use manuscript-review).
license: MIT
metadata: {"version": "1.1.0", "category": "文献与基金", "license": "MIT"}
---

# Deep Critical Paper Reading & Appraisal (1–3 Papers)

Dissect 1 to 3 scientific papers from first principles—reconstructing the **exact claim-to-evidence chain**, testing methodological assumptions, and separating demonstrated findings from author spin.

## 1. Scope & Boundary

- **Use for**: Deep critical reading of **1 paper** (full-text Reading Card + figure-by-figure evidence chain) or head-to-head comparison of **2–3 closely related papers**.
- **Not for (Handoff)**:
  - Cross-study synthesis or evidence matrices across **`>= 5` papers** -> `literature-review`.
  - Pre-submission mock peer review of the user's own unpublished manuscript -> `manuscript-review`.
  - Converting a verified paper into a presentation deck -> `paper-to-slides`.

## 2. Epistemic Contract & Stop Gates

1. **Abstract-Only Degradation Gate (Hard Stop)**:
   - Never fabricate a Fig. 1–5 evidence map, sample sizes, statistical tests, or Methods critique from a title/abstract alone.
   - If full text is NOT provided and cannot be fetched via active tools, **STOP full Reading Card generation** and output strictly an **Abstract-Only Triage Card** with `[FULL TEXT REQUIRED]` markers.
2. **Zero Parametric Hallucination**: Quote or paraphrase only what is present in the retrieved/provided full text, tables, or PDF.
3. **Preprint & Retraction Caution**: Label bioRxiv/medRxiv/arXiv papers `[PREPRINT — NOT PEER-REVIEWED]` unless a published version is verified. Check for errata/retractions when tools allow.

## 3. Phi Runtime & Context Routing

Follow this strict retrieval priority chain before reading:
1. **`literature-library` plugin (if enabled)**: `lit_get(with_fulltext=True)` for known IDs, `lit_search` to locate, `lit_curate(action="versions")` to check preprint-to-published status.
2. **Database MCPs (if enabled)**: `pubmed` / `biorxiv` / `europepmc`.
3. **Local Workspace Files**: User-provided PDF, Markdown, LaTeX, or text files.
4. **Stop Gate Fallback**: If 1–3 do not yield full text -> degrade immediately to **Abstract-Only Triage Card** and request the PDF.

## 4. Core Workflow & Mode Routing

| Input Available | Mode | Allowed Output |
| --- | --- | --- |
| Full text (1 paper) | **Full Critical Reading Card** | Figure-by-figure evidence chain, Methods/stats audit, overclaim check, reproducibility |
| Full text (2–3 papers) | **Comparative Reading Matrix** | Side-by-side design, cohort `n`, metrics, conflicting findings, actionable transfer |
| Title / Abstract only | **Abstract-Only Triage Card** | Claimed advance + 3–5 verification questions to check once PDF is uploaded (`[FULL TEXT REQUIRED]`) |

1. **Verify Source & Status**: Declare Full-Text vs Abstract-Only mode, journal/year, DOI/PMID, and peer-reviewed vs preprint status.
2. **Reconstruct Evidence Chain (Full Text Only)**:
   - Classify study design on the 6-level causal ladder (Descriptive -> Perturbation + Rescue).
   - Map each key figure: `Figure/Panel | Core Claim | System & Biological n | Key Metric (Effect Size, 95% CI, P/FDR) | Method/Control Caveat`.
3. **Stress-Test & Separate Spin from Evidence**: Compare Abstract/Title claims against actual Results, controls, pseudoreplication, external validation, and data/code availability (`GEO`/`Zenodo`/`GitHub`).
4. **Actionable Transfer**: Extract 2–3 reusable methods, design lessons, or testable hypotheses for the user's own work.

## 5. Worked Micro-Example (Abstract-Only Stop Gate)

- *Input*: User pastes only an abstract (`AUC = 0.85` 12-gene signature) and asks for a Fig. 1–4 breakdown.
- *Correct Output*:
  > **[FULL TEXT REQUIRED — Abstract-Only Triage Card]**
  > Full text is unavailable, so Fig. 1–4 and Methods cannot be reconstructed without fabrication.
  > - **Claimed Advance**: 12-gene signature predicting anti-PD-1 response (`AUC = 0.85`).
  > - **Check Once PDF Is Uploaded**: (1) Is `AUC = 0.85` internal CV or independent external validation (`TRIPOD+AI` leakage check)? (2) Exact biological `n` and event count? (3) Are feature-selection scripts and raw data deposited?

## 6. Output Contract

- **Full Text**: (1) Citation & Status Header, (2) 3–4 sentence Executive Summary, (3) Figure-by-Figure Evidence Chain Table, (4) Critical Appraisal (Strengths, Overclaims/Limitations, Reproducibility), (5) Actionable Takeaways (or Comparison Matrix for 2–3 papers).
- **Abstract Only**: `[FULL TEXT REQUIRED]` Abstract-Only Triage Card + targeted verification questions.
