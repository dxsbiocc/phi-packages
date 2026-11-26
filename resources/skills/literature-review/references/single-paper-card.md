# Deep Critical Paper Reading Branch (1–3 Papers)

Load this reference ONLY when reading or comparing 1 to 3 papers in depth.

## 1. Full Critical Reading Card Template (Full Text Required)

```markdown
### Paper Metadata
- **Citation**: Author et al., *Journal* (Year) | DOI/PMID: `...` | Status: `Peer-reviewed` / `[PREPRINT — NOT PEER-REVIEWED]`
- **Design & Causal Level**: `[e.g., Level 2 Observational Cohort + Level 5 In Vitro Perturbation]`

### 1. Core Question & Advance (3–4 sentences)
[Bottleneck -> Approach -> Headline quantitative finding -> Main boundary condition.]

### 2. Figure-by-Figure Evidence Chain
| Figure | Core Claim | System & Biological `n` | Key Quantitative Evidence | Method/Control Caveat |
| --- | --- | --- | --- | --- |
| Fig. 1 | ... | `n = ...` | Effect size, `95% CI`, `P/FDR` | ... |

### 3. Critical Appraisal
- **Strongest Evidence**: ...
- **Overclaim / Causal Gap**: Separate Title/Abstract spin from what Results actually prove.
- **Statistical & Control Audit**: Pseudoreplication, confounders, external validation, negative/positive controls.
- **Reproducibility**: Data `[GEO/Zenodo ID or Missing]` | Code `[GitHub/DOI or Missing]`

### 4. Actionable Transfer / Side-by-Side Comparison (if 2–3 papers)
- ...
```

## 2. Worked Micro-Example (Abstract-Only Stop Gate)

- *Input*: User pastes only an abstract (`AUC = 0.85` 12-gene signature) and asks for a Fig. 1–4 breakdown.
- *Correct Output*:
  > **[FULL TEXT REQUIRED — Abstract-Only Triage Card]**
  > Full text is unavailable via workspace files or retrieval tools, so Fig. 1–4 and Methods cannot be reconstructed without fabrication.
  > - **Claimed Advance**: 12-gene signature predicting anti-PD-1 response (`AUC = 0.85`).
  > - **Verify Once PDF Is Uploaded**: (1) Internal CV vs independent external validation (`TRIPOD+AI` leakage)? (2) Exact biological `n` and event count? (3) Public code & data availability?
