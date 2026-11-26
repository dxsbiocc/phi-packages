---
name: humanizer
description: Rewrite AI-sounding text so it reads like natural human writing without changing meaning, facts, or scientific rigor. Enforces an Academic Register Lock for scientific prose (preserving hedges, causal level, numbers, and LaTeX/citation anchors) while removing AI tells such as not-X-but-Y contrasts, forced triads, inflated hype, em-dash overuse, and staged openers. Use for 去除AI味, 去机器味, 人性化改写, and de-AI prose editing. Not for drafting from raw data (use manuscript-writing) or full academic restructuring/translation/word-count compression (use manuscript-polishing).
license: MIT
metadata: {"version": "1.1.0", "category": "论文写作", "license": "MIT"}
---

# Humanizer (De-AI Prose Editing with Academic Register Lock)

Strip recognizable LLM writing tics so prose sounds like a thoughtful human author—while **never degrading scientific precision, epistemic hedges, or factual accuracy**.

## 1. Scope & Boundary

- **Use for**: Removing LLM stylistic fingerprints from academic manuscripts, grant narratives, cover letters, blog posts, or technical docs.
- **Not for (Handoff)**:
  - Drafting sections from tables/figures or rebuilding manuscript arguments -> `manuscript-writing`.
  - Academic CN <-> EN translation, journal formatting, or word-count compression -> `manuscript-polishing`.
  - Reviewing scientific or statistical validity -> `manuscript-review` / `statistical-review`.

## 2. Epistemic Contract & Stop Gates

1. **Academic Register Lock (Non-Negotiable for Scientific Text)**:
   - Never delete legitimate scientific hedges (*"suggests"*, *"is consistent with"*, *"within this retrospective cohort"*).
   - Never upgrade a correlational verb to a punchier causal verb (*"associates with"* must NEVER become *"drives"* or *"triggers"*).
   - Never inject casual idioms, rhetorical questions, or blog tone into formal journal prose.
2. **Zero Fact/Anchor Drift**: Every number, unit, sample size `n`, `P/FDR`, gene/chemical symbol, and markup anchor (`\cite{}`, `\ref{}`, `[@key]`, Markdown links) must survive verbatim.
3. **Meaning-Loss Stop Gate**: If a paragraph is pure empty AI filler with no concrete underlying facts, **do not invent specifics** to make it sound human. Flag `[AUTHOR INPUT NEEDED: replace generic filler with specific mechanism/result]` instead.

## 3. Phi Runtime & Context Routing

- **Register Auto-Detection**: Inspect the input (and surrounding workspace manuscript/proposal files if present) to lock one of three tiers:
  - `Formal Academic` (journal paper, rebuttal, NSFC/NIH grant) -> **Academic Register Lock ON**.
  - `Semi-Formal Technical` (README, lab report, talk script) -> natural voice, preserve technical specs.
  - `General / Editorial` (blog, outreach, essay) -> fuller rhythmic variation allowed.
- **Standalone**: Zero external tools or sibling skills required.

## 4. Core Workflow & AI-Tell Replacement Table

1. Classify register and lock all numbers, citations, and epistemic hedges.
2. Audit and rewrite concrete AI tells using the table below:

| AI Tell | Typical Pattern (EN / ZH) | Humanized Fix (Academic-Safe) |
| --- | --- | --- |
| **Staged Openers** | *"In recent years..."*, *"近年来，随着……的飞速发展"* | Open directly with the specific biological/technical problem |
| **Inflated Lexicon** | *delve, tapestry, pivotal, nuanced, underscores, multifaceted* / *至关重要, 赋能, 深刻揭示* | Plain, exact verbs (*examine, show, support, remains unclear*) |
| **Not-X-but-Y Crutch** | *"It is not merely X; rather, it is Y."* / *"这不仅是……，更是……"* | State `Y` directly, or state both `X` and `Y` plainly |
| **Forced Triads** | *"robust, scalable, and interpretable"* / *"系统性、创新性和前瞻性"* | Keep only attributes actually evidenced in the text |
| **Em-Dash / Label Spam** | Multiple `—` per paragraph; bold lead-ins inside narrative prose | Integrate into normal sentences with commas or periods |
| **Theatrical Closers** | *"Ultimately, this paves the way for..."* / *"开辟了崭新的道路"* | End on the concrete implication or explicit scope boundary |

3. Vary sentence length and openings naturally; dissolve mechanical `"Furthermore / Moreover / Additionally"` chains.
4. Diff every number, citation key, and causal verb against the source before delivery.

## 5. Worked Micro-Examples

- **English Scientific Discussion (`Formal Academic` Lock ON)**:
  - *Before*: *"In recent years, spatial transcriptomics has emerged as a pivotal tool. Importantly, our findings underscore that fibroblasts are not merely passive bystanders; rather, they play a multifaceted, dynamic, and intricate role in shaping immunity (`r = 0.54`, `P = 0.002`, Fig. 3B). Ultimately, this paves the way for transformative oncology."*
  - *After*: *"Our spatial analysis shows that fibroblast density correlates with T-cell exclusion across the tumor margin (`r = 0.54`, `P = 0.002`, Fig. 3B), suggesting that stromal remodeling is linked to local immune suppression rather than occurring independently. Prospective perturbation studies are needed to test whether targeting fibroblasts restores T-cell infiltration."*
- **Chinese Academic / NSFC Narrative (`Formal Academic` Lock ON)**:
  - *Before*: *"近年来，随着单细胞测序技术的飞速发展，肿瘤研究迎来了前所未有的机遇。本研究不仅系统性地解析了巨噬细胞异质性，更深刻揭示了其在耐药中的关键作用，具有重要的理论意义、临床价值和转化前景。"*
  - *After*: *"单细胞测序使在细胞亚群层面追踪肿瘤耐药成为可能。前期数据表明，耐药样本中 *SPP1*+ 巨噬细胞比例升高（`n = 18` 例），但该亚群是否直接参与耐药维持尚缺乏功能证据。本项目拟通过队列验证与功能敲降检验其作为干预靶点的可行性。"*

## 6. Output Contract

1. **Rewritten Text** (matching input format and register).
2. **Register Applied** (`Formal Academic [Register Lock ON]` / `Semi-Formal Technical` / `General`).
3. **Brief Edit Note** (2–4 bullets listing removed AI patterns and confirming zero drift in facts, causal level, or citations).
