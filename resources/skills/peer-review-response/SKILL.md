---
name: peer-review-response
description: Draft strategic, evidence-backed point-by-point responses to journal peer reviewers, triage comments by feasibility, map every change to exact manuscript locations, and write the resubmission cover letter. Use for 审稿回复, 修回信, rebuttal letter, response to reviewers, and major/minor revision planning. Not for pre-submission mock review (use manuscript-review), initial-submission cover letters (use manuscript-writing), or standalone statistical audits (use statistical-review).
license: MIT
metadata: {"version": "1.1.0", "category": "论文写作", "license": "MIT"}
---

# Peer-Review Response & Rebuttal Strategy

Turn critical peer-review reports into a calm, structured, evidence-backed **Point-by-Point Response Letter** and **Resubmission Cover Letter** that make every revision effortless for editors and reviewers to verify.

## 1. Scope & Boundary

- **Use for**: Triaging received reviewer/editor comments, drafting verbatim point-by-point responses with exact manuscript locations and pasted excerpts, and writing the **Resubmission Cover Letter**.
- **Not for (Handoff)**:
  - Pre-submission mock peer review -> `manuscript-review`.
  - Initial-submission cover letter -> `manuscript-writing`.
  - Independent statistical audit of raw data -> `statistical-review`.
  - Full manuscript polishing -> `manuscript-polishing`.

## 2. Epistemic Contract & Stop Gates

1. **Zero Phantom Revisions**: Never claim *"We performed the new analysis/experiment and found `P = 0.01`"* unless those results exist in verified user files or input. Mark uncompleted work `[PENDING ACTION: <task>]` and unknown anchors `[LINE/PAGE TBD]`.
2. **Verbatim Reviewer Integrity**: Split compound paragraphs into numbered atomic items (`R1.1`, `R1.2`, ...), quoting every sub-question verbatim—never quietly drop or merge away a hard question.
3. **Unresolved Data Stop Gate**: If reviewers demand new analyses/experiments not yet run, deliver the **Triage Matrix** + draft response skeletons with explicit `[PENDING RESULT: ...]` slots rather than faking completed results.

## 3. Phi Runtime & Context Routing

- **Workspace Diff Grounding**: When original/revised manuscripts, new tables, or `*.phi-artifact.json` sidecars are in the workspace, inspect them to cite exact sections, figure panels, and verified numbers.
- **Standalone**: Works with zero external plugins or sibling skills.

## 4. Core Workflow & Triage Strategy Matrix

| Comment Type | Strategy | Mandatory Evidence |
| --- | --- | --- |
| **Valid Gap — Feasible** | Accept & Execute | Exact new metric (`n`, effect size, CI, `P/FDR`) + new panel + quoted text |
| **Valid Critique — Out of Scope / Infeasible** | Acknowledge + Scope Boundary + Partial Remedy | Candid Limitation + calibrated claim strength + best feasible proxy/sensitivity check |
| **Reviewer Misread Existing Text** | Clarify Without Blame | Take responsibility for visibility + quote existing data + add main-text signposting |
| **Methodologically Flawed Request** | Polite Technical Pushback | Concise statistical/biological rationale + why current method is valid |
| **Minor / Editorial** | Accept Directly | 1-sentence confirmation + exact location |

For every atomic comment, use the **4-Part Response Anatomy**:
1. **Verbatim Comment** (`R#.#`)
2. **Direct Answer** (1–2 sentences: what was done or why you respectfully disagree)
3. **Evidence / Rationale** (verified numbers or methodological reasoning)
4. **Changes Made** (`Section, Page, Lines, Figure` + exact pasted revised text in block quotes)

Tone rules: Never defensive (*"The reviewer misunderstood..."*) and never sycophantic (no *"brilliant insight"* repeated 20 times); keep revisions mutually consistent across Reviewer 1/2/3.

## 5. Worked Micro-Example (Polite Pushback + Claim Calibration)

- **R2.4**: *"To prove Gene X causes resistance, generate a conditional KO mouse model."*
- **Response**: *"We agree that in vivo perturbation is required to establish causal necessity. Generating a conditional KO line is beyond the scope of this clinical-observational study. To keep conclusions strictly calibrated to our evidence, we have (1) revised all causal verbs to associative/prognostic phrasing, (2) added an explicit Limitations paragraph, and (3) added an independent external cohort (`n = 86`, `AUC = 0.81`, `95% CI [0.72, 0.90]`, Fig. 4D)."*
- **Changes made** (`Discussion, p. 16, lines 382–388; Fig. 4D`): *"> Because our cohorts are observational, these findings establish Gene X as a reproducible clinical correlate of resistance rather than a proven causal driver..."*

## 6. Output Contract

1. **Revision Triage Table** (`ID | Severity | Strategy | Status [Done/Pending] | Location`).
2. **Resubmission Cover Letter** (3–5 thematic major revisions for the editor + neutral note on any out-of-scope items).
3. **Point-by-Point Response Letter** (every atomic comment with 4-part anatomy).
4. **Pending Actions & Anchor Manifest** (`[PENDING ACTION]` and `[LINE/PAGE TBD]` items).
