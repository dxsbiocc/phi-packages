---
name: manuscript-polishing
description: Polish, tighten to word limits, translate (Chinese-English), or clean up LaTeX/Markdown formatting of existing academic manuscripts while strictly preserving facts, causal strength, technical terminology, and citations. Use for 论文润色, 学术翻译, 正文精简, 字数压缩, sentence flow, and journal style alignment. Not for drafting new sections from raw data (use manuscript-writing), pure AI-tell de-homogenization (use humanizer), mock peer review (use manuscript-review), or rebuttal drafting (use peer-review-response).
license: MIT
metadata: {"version": "1.1.0", "category": "论文写作", "license": "MIT"}
---

# Academic Manuscript Polishing, Trimming & Translation

Improve clarity, conciseness, and idiomatic academic prose of existing manuscripts **without altering factual content, statistical precision, or causal claim strength**.

## 1. Scope & Boundary

- **Use for**: Line-level academic polishing (`Light` / `Standard` / `Structural`), disciplined word-count compression (`Word-Count Trim`), CN <-> EN academic translation, and LaTeX/Markdown typography cleanup.
- **Not for (Handoff)**:
  - Drafting new sections from raw data or rebuilding the whole argument -> `manuscript-writing`.
  - Removing AI prose tells as the primary goal -> `humanizer`.
  - Pre-submission mock review -> `manuscript-review`.
  - Post-decision rebuttal -> `peer-review-response`.

## 2. Epistemic Contract & Stop Gates

1. **Fact & Number Invariance**: Never alter any number, unit, `P/FDR`, CI, sample size `n`, gene/protein symbol, reagent, or citation key.
2. **Causal Strength Lock**: Never upgrade *"associated with"* to *"drove"* or *"caused"* for stylistic punch.
3. **Ambiguity Stop Gate**: If a source sentence is scientifically ambiguous (e.g., unclear whether `n = 12` refers to mice or slices), **do NOT guess**—insert `[AUTHOR CHECK: <specific ambiguity>]` inline and list it in the Change Log.

## 3. Phi Runtime & Context Routing

- **Workspace Terminology Lock**: Inspect existing manuscript files, `.bib`, or `.phi-artifact.json` (if present) to preserve exact gene casing, abbreviations, and US vs UK spelling.
- **Exact Word Counting**: When tightening to a hard word limit (`<= 150` or `<= 200` words), verify word count via Python/shell if available.
- **On-Demand Branch References** (load ONLY if the user's mode matches):
  - **CN <-> EN Academic Translation** -> read [references/zh-en-translation.md](references/zh-en-translation.md).
  - **LaTeX / Markdown / Biomedical Typography Cleanup** -> read [references/latex-typography.md](references/latex-typography.md).

## 4. Core Workflow & Mode Matrix

| Mode | Scope of Edits | What Stays Untouched |
| --- | --- | --- |
| **Light Copyedit** | Grammar, articles, prepositions, punctuation, tense consistency | Sentence order, paragraph structure, author voice |
| **Standard Academic** | Clarity, concision, active verbs, given-before-new cohesion | Paragraph logic, all data/citations, causal level |
| **Structural Flow** | Reordering sentences within paragraphs, sharpening topic sentences | Core scientific claims, evidence scope |
| **Word-Count Trim** | 15%–35% compression to hit journal limits (trim fluff, never metrics) | `n`, effect size, CI, `P/FDR`, key limitations |
| **Academic Translation** | Idiomatic CN <-> EN translation (load `references/zh-en-translation.md`) | Exact values, citation keys, causal scope |

1. Confirm mode, spelling variant (US/UK), and any hard word limit.
2. Lock all numbers, symbols, and markup anchors (`\cite{}`, `\ref{}`, `[@key]`).
3. Run 3 passes: (1) Precision & Causal Lock, (2) Concision & Active Flow, (3) Tense & Format (Past for your Methods/Results; Present for established facts and figure references).
4. Diff all numbers, units, and citation keys against source before delivery.

## 5. Worked Micro-Example (Word-Count Trim)

- *Before (66w)*: *"In this study, we developed a deep learning framework named CellMap-Net for the purpose of predicting spatial gene expression directly from routine H&E histology whole-slide images. Across three independent external cohorts comprising a total of 412 patients, our proposed method achieved a mean Pearson correlation coefficient of 0.68 (`SD = 0.07`), outperforming baseline models by 14.2% (`P < 0.001`)."*
- *After (43w — zero metric loss)*: *"We developed CellMap-Net, a deep learning model that predicts spatial gene expression from H&E whole-slide images. Across three external cohorts (`n = 412` patients), CellMap-Net achieved mean Pearson `r = 0.68` (`SD = 0.07`), outperforming baselines by 14.2% (`P < 0.001`)."*

## 6. Output Contract

1. **Polished / Translated Text** (drop-in ready, matching input format).
2. **Word Count Check** (`Before -> After`, if a limit was requested).
3. **Editorial Change Log** (3–6 concise bullets + any `[AUTHOR CHECK]` items).
