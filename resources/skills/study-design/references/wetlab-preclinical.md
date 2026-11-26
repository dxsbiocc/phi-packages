# Wet-Lab & Preclinical Animal Experimental Design Branch

Load this reference ONLY when designing cell-culture, biochemical, organoid, or preclinical animal experiments.

## 1. Control & Experimental Unit Matrix

| Assay / Model | True Biological Unit (`n`) | Common Pseudoreplication Trap | Mandatory Controls |
| --- | --- | --- | --- |
| **Cell culture** | Independent biological passage/donor line or independent differentiation batch | Treating 3 technical wells on the same plate as `n = 3` | Vehicle/DMSO, non-targeting sgRNA/shRNA, positive control, rescue construct |
| **Organoids / PDOs** | Independent patient donors (`n` donors) | Pooling 50 organoids from 1 donor as `n = 50` | Matched donor vehicle + `>= 3` independent donors |
| **Animal (individual dosing)** | Individual animal (stratified by litter & sex) | Measuring 10 histological sections per mouse as `n = 10` | Vehicle/sham, littermate wild-type/Cre-only control, blinded scoring |
| **Animal (cage dosing / microbiome)** | **Cage** (co-housed mice share microbiota & water/diet intake) | Treating 5 mice in 1 cage as `n = 5` independent units | Multiple cages per arm (`>= 3–5` cages/arm); co-housing or fecal transplant control |

## 2. ARRIVE 2.0 Pre-Registration Essentials

- **Randomization**: Block-randomize by **litter, baseline weight, and sex** using a computer-generated seed (never grab the first 6 mice from a cage).
- **Blinding**: Separate coded drug preparation, outcome measurement/scoring, and primary statistical analysis.
- **Exclusion Criteria**: Pre-specify humane/technical exclusions (e.g., failed injection logged *before* outcome readout).

## 3. Worked Micro-Example (Strong Inference + De-Confounded Design)

- *Strong Inference Split*: Bulk/scRNA-seq shows *SPP1*+ macrophages correlate with T-cell exclusion (`r = 0.61`). `H0` = hypoxia passenger; `H1` = direct SPP1-CD44 T-cell inhibition; `H2` = paracrine myofibroblast collagen barrier.
- *Crucial Experiment*: 2x2 co-culture (CD8+ T cells +/- *SPP1*-KO macrophages +/- myofibroblasts) under matched normoxia/hypoxia across `n = 6` independent donor-derived lines (balanced across plate columns and assay days). `H1` predicts rescue without myofibroblasts; `H2` predicts rescue ONLY with myofibroblasts; `H0` predicts null effect under matched oxygen.
