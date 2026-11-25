# Mendelian Randomization & Colocalization Playbook (STROBE-MR)

## 1. GWAS QC & Clumping (`mr_gwas_qc`)

- **Check $\lambda_{\text{GC}}$**: $0.95 \le \lambda_{\text{GC}} \le 1.10$ is typical for well-controlled cohorts; $\lambda_{\text{GC}} > 1.20$ without polygenic LDSC intercept confirmation warns of population stratification or cryptic relatedness.
- **Clumping**: when an external 1000G PLINK reference panel is unavailable, `clump_kb=10000` (10 Mb physical window) greedily retains the lowest-$P$ variant per window to suppress cis-LD correlation; state this explicitly in methods.

## 2. Allele Harmonization & Instrument Strength (`mr_harmonize`)

- **Weak-instrument red line**: require per-SNP $F = \beta_{\text{exp}}^2 / \text{SE}_{\text{exp}}^2 \ge 10$.
- **Palindromic SNPs (`A/T` and `C/G`)**: forward vs reverse strand is indistinguishable near $\text{EAF} = 0.50$. Drop any palindromic variant with $0.42 \le \text{EAF} \le 0.58$ (or missing EAF).

## 3. Evaluating Causal Evidence (`mr_estimate`)

- **Primary estimator**: Multiplicative random-effects IVW ($K \ge 2$) or Wald ratio ($K = 1$).
- **Triangulation across robust estimators ($K \ge 3$)**:
  1. Consistent effect sign across **IVW**, **Weighted Median**, and **MR-Egger**;
  2. **MR-Egger intercept $p \ge 0.05$** (no evidence of directional horizontal pleiotropy);
  3. **Leave-One-Out (`mr_singlesnp_loo.tsv`)**: no single pleiotropic locus (e.g. *APOE* or *HLA*) drives the entire IVW significance;
  4. **Cochran's $Q$**: when $p_Q < 0.05$, inspect `outlier_snps` and rely on multiplicative random-effects IVW + Weighted Median rather than fixed-effects IVW.

## 4. Cis-MR & Bayesian Colocalization (`mr_coloc`)

- In drug-target cis-pQTL / cis-eQTL MR ($K = 1\text{–}3$ variants), a significant Wald ratio frequently arises from LD between two distinct causal variants (`H3`).
- Interpret posterior probabilities on the $\pm 250\text{–}500\text{ kb}$ cis window (`p1=1e-4, p2=1e-4, p12=1e-5`):
  - `PP.H4 >= 0.80`: strong evidence of a shared causal variant.
  - `PP.H3 >= 0.80`: strong evidence of two distinct LD-linked variants (confounding by LD — **reject** shared causal mechanism!).
