# Journal Tier Presets & Initial Cover Letter Guide

Load this reference when targeting a specific journal format (*Nature*, *Science*, *Cell*, clinical journals, or bioinformatics/methods journals) or drafting an initial-submission Cover Letter.

## 1. Journal Tier & Structural Matrix

| Target Tier | Abstract Limit | Intro Length | Results & Discussion | Key Format Expectations |
| --- | --- | --- | --- | --- |
| *Nature* Article | `<= 200`w unreferenced summary paragraph | 3–4 short paragraphs (~500w) | Separate Results & Discussion | `<= 4–6` main figures; online Methods; Reporting Summary |
| *Science* Research Article | `<= 125`w | 2–3 tight paragraphs | Tightly linked or combined | `<= 4–5` main figures; `~2,500–4,500`w main text |
| *Cell* Article | `<= 150`w + 4 Highlights + eTOC blurb | 3–5 paragraphs | Separate Results, Discussion, Limitations of the Study | STAR Methods + Key Resources Table (RRIDs mandatory) |
| *Nature Communications* / *Science Advances* | `<= 150`w | 3–4 paragraphs | Separate Results & Discussion | Flexible main figures (`<= 8`); explicit Data/Code Availability |
| Clinical (*Lancet*, *NEJM*, *JAMA*, *BMJ*) | `<= 250–300`w structured (Background/Methods/Findings/Interpretation) | 2–3 crisp paragraphs | Strict Results (no interpretation) + structured Discussion | Trial registration, primary/secondary endpoints, absolute + relative risk |
| Bioinformatics / Methods (*Nature Methods*, *Genome Biology*, *Bioinformatics*) | `<= 150–250`w | 3–4 paragraphs | Benchmark + Ablation + Real-data Discovery | Runtime/memory scaling, public reproducible repo, comparison vs SOTA baselines |

## 2. Worked Before/After Micro-Examples

- **Results Subsection (Causal Calibration + Donor-Level `n`)**:
  - *Before*: *"Strikingly, *SPP1*+ macrophages drive fibrotic remodeling in non-responders (`P < 1e-300`, Wilcoxon across 42,000 cells, Fig. 2A)."*
  - *After*: *"*SPP1*+ macrophages were enriched in non-responders relative to responders (median 24.3% vs 7.8% of myeloid cells across `n = 18` vs `n = 14` donors; pseudobulk `log2FC = 1.64`, `95% CI [0.92, 2.36]`, `FDR = 0.004`, Fig. 2A–B) and correlated with stromal fibrosis score (`Spearman rho = 0.61`, `n = 32` donors, `P = 2.1e-4`, Fig. 2C)."*
- **Discussion Opening (Advance vs Repetition)**:
  - *Before*: *"In this study, we found Gene X was upregulated 2.3-fold (`P = 0.003`). Future larger cohorts are needed."*
  - *After*: *"Our findings position Gene X as a reproducible transcriptional correlate of early relapse across two independent cohorts, reconciling prior bulk-tissue discrepancies `[CITE: prior bulk study]`. Because our design is observational, CRISPRi knockdown in patient-derived organoids will directly test causal necessity."*

## 3. Initial-Submission Cover Letter Blueprint

1. **Opening (3–4 sentences)**: Title, target journal, one-sentence concrete advance + headline quantitative result.
2. **Why Now & Why This Journal (1 paragraph)**: Exact bottleneck solved and why the journal's readership cares.
3. **Integrity & Logistics (bullets)**: Originality confirmation, preprint DOI (if any), competing interests, Data/Code repository status, optional suggested/excluded reviewers.
