# Multi-Paper Systematic & Scoping Synthesis Branch (>= 5 Studies)

Load this reference ONLY when synthesizing `>= 5` studies into a thematic review and Evidence Matrix.

## 1. Evidence Matrix Schema

| Study (Author, Year, ID) | Design & Causal Level | System / Cohort (`n`) | Primary Metric (Effect Size, CI, `P`) | Key Limitation / Risk of Bias | Status |
| --- | --- | --- | --- | --- | --- |
| Chen 2023 (`PMID:...`) | Retrospective cohort (L2–3) | `n = 142` CRC (2 centers) | `HR = 1.78`, `95% CI [1.21, 2.62]` | No external validation; bulk RNA | Peer-reviewed |
| Li 2024 (`bioRxiv:...`) | scRNA-seq + mouse KO (L5) | `n = 18` donors; `n = 6`/arm | `log2FC = 1.52`, `FDR = 0.008` | Small mouse arm; unblinded | `[PREPRINT]` |

## 2. Thematic Synthesis Rules

- **Never write chronological laundry lists** (*"Smith found X. Jones found Y. Wang found Z."*). Group by **mechanism, methodology, or clinical axis**.
- **Reconcile Contradictions Explicitly**: Explain why studies disagree (disease stage, assay platform, sample size/event power, confounding adjustment, biological unit).
- **Weight Evidence by Hierarchy**: Multi-center RCT / prospective / external-validation studies outweigh small single-center or cell-line-only studies; note publication bias and negative-result gaps.

## 3. Worked Micro-Example (Thematic Reconciliation)

- *Strong Synthesis*: *"Across four GI tumor cohorts (`k = 4`, `n = 518`), elevated *SPP1*+ macrophage abundance associated with shorter OS in three studies (`HR` `1.62–2.14`; Smith 2021 `PMID:...`, Jones 2022 `PMID:...`, Park 2023 `PMID:...`), whereas Wang 2023 (`PMID:...`, `n = 46`) reported a null result (`HR = 1.12`, `95% CI [0.71, 1.77]`), likely reflecting Stage I–II enrollment and low event power (`19` deaths). Secondary review Brown 2024 was excluded from quantitative pooling."*
