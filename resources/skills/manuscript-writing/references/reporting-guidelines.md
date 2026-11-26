# Study-Design Reporting Guidelines & Integrity Checklist

Load this reference when drafting or auditing Methods, Results compliance, Ethics, or Data/Code Availability for a specific study design.

## 1. Mandatory Reporting Guideline Routing Table

| Study Design | Standard | Mandatory Items to Include in Methods & Results |
| --- | --- | --- |
| **Randomized Controlled Trial** | **CONSORT 2025** | Trial registration ID, randomization sequence, allocation concealment, blinding, ITT vs per-protocol, participant flow, harms |
| **Observational (Cohort / Case-Control / Cross-Sectional)** | **STROBE** | Eligibility criteria, confounder selection rationale, missing-data handling, sensitivity analyses, follow-up completeness |
| **Animal / Preclinical** | **ARRIVE 2.0** | Species/strain/sex/age, randomization, outcome-assessor blinding, sample-size justification, exclusions, welfare/IACUC |
| **Clinical Prediction / ML / AI** | **TRIPOD+AI** | Leak-free train/val/test split, external validation cohort, discrimination (`AUC`/`PR-AUC`) + calibration, class balance, code availability |
| **Diagnostic Accuracy** | **STARD** | Reference standard, pre-specified cutoff, sensitivity/specificity with `95% CI`, STARD flow |
| **Systematic Review / Meta-Analysis** | **PRISMA 2020** | Databases, exact query date & string, screening counts, risk-of-bias tool, heterogeneity (`I^2`, `tau^2`) |
| **Mendelian Randomization** | **STROBE-MR** | Instrument selection (`F > 10`), pleiotropy sensitivity tests (MR-Egger, weighted median, MR-PRESSO), sample overlap check |

## 2. Universal Biomedical & Computational Integrity Items

- **SAGER (Sex & Gender)**: Report sex/gender composition, disaggregate primary outcomes by sex where applicable, or justify single-sex design.
- **Reagents & Cell Lines**: STR profiling and mycoplasma testing dates/status for cell lines; antibody/reagent catalog numbers and RRIDs.
- **Image Integrity**: Full uncropped Western blots (with molecular weight markers) in Source Data; microscopy scale bars and uniform linear adjustments only.
- **AI Disclosure**: Disclose any generative AI assistance per journal policy without attributing authorship.
- **Data & Code Availability**: Persistent accessions (`GEO`, `EGA`, `PRIDE`, `Zenodo` DOI) and pinned software/environment versions (`R`, `Python`, key packages).
