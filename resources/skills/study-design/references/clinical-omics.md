# Clinical Cohort & High-Throughput Omics Design Branch

Load this reference ONLY when designing human clinical studies, biomarker discovery/validation cohorts, or bulk/single-cell/spatial omics studies.

## 1. Omics & Clinical De-Confounding Matrix

| Study Type | Primary Confounding Risk | Mandatory Design Fix |
| --- | --- | --- |
| **Bulk / scRNA-seq / Spatial** | All cases on Run 1 / Chip 1, all controls on Run 2 / Chip 2; unequal dissociation times | **Balanced multiplexing** (hash-tagging/genotype demux) or stratified allocation so every chip/run contains balanced cases, controls, and sexes |
| **Biomarker / ML Discovery** | Splitting slides/samples from the same patient across train & test; overfitting with `EPV < 10` | Pre-allocate an **independent external validation center/cohort**; enforce patient-level splits and `>= 10–20` events per candidate predictor |
| **Survival Cohort** | Short follow-up yielding too few events (power depends on **event count**, not total `N`) | Size the cohort by **required events** (Schoenfeld formula) + expected accrual/follow-up duration and 15% loss-to-follow-up buffer |
| **Longitudinal / Paired** | Dropouts over time; carryover effects in crossover designs | Paired donor design (removes inter-donor variance); plan washout period + mixed-effects model (`LMM`) |

## 2. Single-Cell Power Reality Check

- Increasing cells per donor from `2,000` to `20,000` improves rare-cell detection, **NOT** donor-level differential expression power.
- For case-control pseudobulk DE, prioritize **more biological donors** (`>= 6–10` donors/arm for strong effects; `>= 15–30`/arm for heterogeneous human cohorts) over ultra-deep sequencing of 2 donors.

## 3. Worked Micro-Example (Multicenter Biomarker Cohort)

- *Goal*: Validate a binary treatment-response signature (`expected AUC = 0.80`, `35%` response rate, 3 centers).
- *Design*: Reserve Center C (`n = 80`, ~28 responders) untouched as **external validation**; lock preprocessing and cutoff on Centers A+B (`n = 160`) via nested patient-level CV; balance RNA extraction batches across centers and response status.
