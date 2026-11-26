# Chinese-English Academic Translation Reference

Load this reference ONLY when performing `CN -> EN` or `EN -> CN` academic translation.

## 1. Common CN -> EN False Friends & Chinglish Fixes

| Literal CN->EN Pattern | Problem | Preferred Academic English |
| --- | --- | --- |
| *"With the rapid development of..."* / *"In recent years..."* | Empty throat-clearing | Open directly with the specific scientific gap |
| *"obviously / dramatically / significantly increased"* (no stats) | *"significantly"* is reserved for statistical tests | *"markedly increased"* or *"significantly increased (`P = 0.002`)"* |
| *"prove / verify / confirm"* (from observational data) | Overstates causal certainty | *"support"*, *"indicate"*, *"are consistent with"* |
| *"carry out / conduct / perform an analysis of..."* | Nominalization bloat | *"analyzed"*, *"quantified"*, *"evaluated"* |
| *"the results showed that..."* repeated every sentence | Weak grammatical subject | Make the biological/methodological finding the subject |
| *"has great potential to..."* | Vague hype | *"may enable..."* or *"provides a candidate target for..."* |

## 2. Worked LaTeX-Safe CN -> EN Micro-Example

- **Source (Chinese + LaTeX)**:
  > 如图~\ref{fig:surv}所示，在校正年龄、性别和肿瘤分期后，高表达组的总生存期显著短于低表达组（校正后 `HR = 1.82`，`95% CI 1.24–2.67`，`P = 0.002`，详见文献~\cite{smith2023}的队列定义）。
- **Target (Idiomatic Academic EN, anchors intact)**:
  > As shown in Fig.~\ref{fig:surv}, after adjusting for age, sex, and tumor stage, the high-expression group had significantly shorter overall survival than the low-expression group (adjusted `HR = 1.82`, `95% CI [1.24, 2.67]`, `P = 0.002`; cohort definitions follow~\cite{smith2023}).

## 3. EN -> CN Academic Rules

- Preserve standard Chinese discipline terminology (e.g., `pseudobulk` -> 伪批量, `confounder` -> 混杂因素, `hazard ratio` -> 风险比, `competing risk` -> 竞争风险).
- Keep gene/protein symbols, statistical symbols (`HR`, `95% CI`, `FDR`), and software names in English; avoid clunky literal translations of long English clauses by breaking into 2 crisp Chinese sentences.
