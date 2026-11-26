---
name: paper-to-slides
description: Convert a verified scientific paper, manuscript, or analysis report into a visual-first presentation deck (group meeting, conference talk, thesis defense, or journal club) with timed speaker notes and Q&A backup slides. Defaults to zero-dependency Marp Markdown or Reveal HTML; emits python-pptx scripts only after pre-checking the Python environment. Use for 论文转PPT, 组会汇报, 学术答辩, 会议报告, presentation deck, and slide outline. Not for drafting manuscripts (use manuscript-writing) or deep critical paper appraisal (use literature-review).
license: MIT
metadata: {"version": "1.1.0", "category": "论文写作", "license": "MIT"}
---

# Scientific Paper & Project to Presentation Slides

Transform verified manuscripts or computational analysis reports into crisp, visual-first scientific presentations without losing quantitative accuracy or depending on unavailable office plugins.

## 1. Scope & Boundary

- **Use for**: Converting a verified paper, manuscript draft, or analysis report into a timed talk deck (5 / 10 / 15 / 30 min) with assertion-evidence titles, layout hints, speaker notes, and Q&A backup slides.
- **Not for (Handoff)**:
  - Drafting or restructuring the manuscript itself -> `manuscript-writing`.
  - Deep critical appraisal of an unfamiliar paper before presenting -> `literature-review` (Deep Reading Mode) first.
  - Statistical validity audit -> `statistical-review`.

## 2. Epistemic Contract & Stop Gates

1. **Zero Fabricated Visuals or Numbers**: Every figure path, `n`, effect size, CI, and `P/FDR` must trace to verified workspace files or user input. Never invent placeholder charts or fake summary numbers.
2. **Missing Source Stop Gate**: If neither source text/PDF nor workspace figures/tables are provided, **STOP** and request source files (or deliver an explicit skeleton marked `[INSERT FIGURE: ...]` and `[TBD: ...]`).
3. **Causal Honesty**: Assertion titles must match the study's exact causal level (never upgrade observational associations into *"drives"* on slides).

## 3. Phi Runtime & Context Routing

- **No Assumed Office Plugin or Managed-Env Mutation**: Packaged Phi apps do not bundle an `office` plugin, and managed base environments (`phi:python@1`) are read-only and may lack `python-pptx` / `marp-cli`.
- **Format Routing**:
  - **Default (Zero-Dependency)**: **Marp Markdown** (`slides.md` with `marp: true` frontmatter and relative workspace figure paths) or **Reveal.js HTML**.
  - **Editable `.pptx` Request**: Pre-check `python -c "import pptx"` first. If available, generate `.pptx`; if unavailable, deliver `slides.md` immediately + standalone `generate_deck.py` with an explicit note (never `pip install` into read-only managed envs).
- **Workspace Scan**: Inspect `*.phi-artifact.json` sidecars and `figures/` / `results/` when present to pull verified figure paths and sample sizes.

## 4. Core Workflow & Pacing Matrix

| Talk Format | Duration | Main Slides | Backup | Structural Allocation |
| --- | --- | --- | --- | --- |
| **Lightning Talk** | 5 min | 5–6 | 2 | 1 Title, 1 Gap, 1 Method, 2 Results, 1 Takeaway |
| **Conference Talk** | 10–12 min | 10–12 | 3 | 1 Title, 2 Gap, 1–2 Design, 4–5 Results, 1–2 Limits, 1 Summary |
| **Group Meeting / JC** | 15–20 min | 14–18 | 3–4 | 1 Title, 2 Context, 3 Methods, 6–8 Panel Results, 2 Critique, 1 Q&A |
| **Thesis Defense** | 30–40 min | 24–30 | 5–8 | 1 Title, 4 Gap, 3 Framework, 12–16 Multi-Aim Results, 3 Synthesis, 1–2 Next |

1. Confirm audience, duration, language, and format (`Marp Markdown` default).
2. Map each slide to **one core assertion (`<= 12` words) + one verified figure/panel** (`<= 35–50` words body, max 3–4 bullets, tables `<= 4x5`).
3. Put detailed explanations and pacing (`~60–90s/slide`) in `<!-- Speaker notes: ... -->`.
4. Append 2–3 data-backed **Q&A Backup Slides** covering assumptions, limitations, and sensitivity checks.

## 5. Worked Marp Slide Micro-Example

```markdown
---
## *SPP1*+ macrophages are enriched in non-responders and correlate with fibrosis

<!-- Layout: Left (60%) figures/fig2_spp1.png [panels A-C] | Right (40%) 3 bullets -->

- **3.1x higher abundance** in non-responders (24.3% vs 7.8%; `n = 18` vs `14` donors)
- **Donor-level pseudobulk**: `log2FC = 1.64`, `95% CI [0.92, 2.36]`, `FDR = 0.004`
- **Stromal link**: correlates with fibrosis (`Spearman rho = 0.61`, `P = 2.1e-4`)

*Source: Fig. 2A–C (observational cohort; perturbation tested in Aim 3).*

<!-- Speaker notes (~75s): Walk Fig. 2A left to right; note donor-level aggregation (n = 32) avoids cell-level pseudoreplication; keep claim associative here. -->
```

## 6. Output Contract

1. **Deck Plan & Figure Inventory** (audience, duration, slide count, verified figure paths).
2. **Complete Slide Deck** (Marp `slides.md` default, layout hints, quantitative footer callouts, speaker notes).
3. **Q&A Backup Slides** (2–3 data-backed slides) + **Missing Asset Manifest** (`[INSERT FIGURE]` / `[TBD]`).
