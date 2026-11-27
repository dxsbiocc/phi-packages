# Statistical Genetics & Post-GWAS Analysis Playbook

## 1. GWAS Summary-Statistics Harmonization (`gwas_sumstats_harmonize`)

- **Supported Caller Schemas**: PLINK2 (`--glm`), REGENIE (`LOG10P`, `ALLELE0`/`ALLELE1`), SAIGE (`Allele1`/`Allele2`, `AF_Allele2`), BOLT-LMM (`P_BOLT_LMM`), and GWAS Catalog harmonized TSVs.
- **Scale Standardization**:
  - Odds ratios are converted via $\beta = \ln(\text{OR})$.
  - $-\log_{10}(P)$ values (`LOG10P`) are converted via $P = 10^{-\text{LOG10P}}$ with $Z$-score underflow recovery $P = 2 \Phi(-|Z|)$.
  - Effective sample size for case-control studies:
    $$N_{\text{eff}} = \frac{4}{\frac{1}{N_{\text{cases}}} + \frac{1}{N_{\text{controls}}}}$$
- **Genomic Inflation**:
  $$\lambda_{\text{GC}} = \frac{\operatorname{median}(\chi^2)}{0.4549364}, \qquad \lambda_{1000,\text{CC}} = 1 + 500(\lambda_{\text{GC}} - 1)\left(\frac{1}{N_{\text{cases}}} + \frac{1}{N_{\text{controls}}}\right)$$
  For quantitative traits without case/control counts, $\lambda_{1000} = 1 + (\lambda_{\text{GC}} - 1)\frac{1000}{N_{\text{median}}}$.
- **Lead Locus Extraction**: Greedy minimum-$P$ distance clumping with minimum separation `lead_window_kb` (default `500 kb`) below `p_threshold` (default `5e-8`).

## 2. LD Score Regression (`gwas_ldsc_interpret`)

- **Univariate LDSC (Bulik-Sullivan et al. 2015)**:
  $$E[\chi^2_j \mid \ell_j] = \frac{N_j h^2}{M}\ell_j + N_j a + 1$$
  Heteroscedasticity weights $w_j = \left[\max(1, \ell_j)\left(1 + \frac{N_j \hat{h}^2 \ell_j}{M}\right)^2\right]^{-1}$ are updated in two stages.
- **Confounding vs. Polygenicity**:
  $$\text{Attenuation Ratio} = \frac{\text{intercept} - 1}{\overline{\chi^2} - 1}$$
  An attenuation ratio near $0$ (with intercept $\approx 1.0$) confirms polygenic signal rather than stratification.
- **Liability-Scale $h^2$ (Lee et al. 2011)**:
  Given sample case proportion $P$ and population prevalence $K$, with threshold $t = \Phi^{-1}(1 - K)$ and density $z = \phi(t)$:
  $$h^2_{\text{liab}} = h^2_{\text{obs}} \cdot \frac{K^2(1 - K)^2}{P(1 - P)z^2}$$
- **Genomic Block Jackknife & Bivariate $r_g$**:
  Contiguous genomic blocks (`n_blocks`, default `50`) yield leave-one-block-out SEs for $h^2$, intercept, and cross-trait genetic correlation $r_g = \frac{\rho_g}{\sqrt{h^2_1 h^2_2}}$.
- **Local LD Proxy Fallback**: When neither an `l2` column nor `--ld_scores` is provided, a local physical-window LD proxy is computed with a warning; interpret absolute $h^2$ cautiously without an ancestry-matched reference panel.

## 3. Fine-Mapping (`gwas_finemap`)

- **SuSiE-RSS (Wang et al. 2020; Zou et al. 2022)**:
  Uses summary $Z$-scores, sample size $N$, and regularized LD matrix $\tilde{R} = (1 - \lambda)R + \lambda I$ (automatically ensuring positive semi-definiteness) via Iterative Bayesian Stepwise Selection (IBSS) with $L \in [1, 10]$ single-effect components.
- **Posterior Inclusion Probabilities & Purity-Filtered Credible Sets**:
  $$\text{PIP}_j = 1 - \prod_{l=1}^L (1 - \alpha_{l,j})$$
  Each component $l$ with positive prior variance $V_l > 0$ forms a minimal `coverage` (default `0.95`) credible set $S_l$ and is retained only if $\min_{i,j \in S_l} |R_{i,j}| \ge \text{min\_abs\_corr}$ (default `0.50`).
- **Single-Causal Fallback**: When `ld_matrix` is omitted, exact Wakefield (2009) Approximate Bayes Factor ($L=1$) fine-mapping is performed with a warning.

## 4. Cross-Cohort Meta-Analysis (`gwas_meta_qc`)

- **Allele Alignment & EAF QC**: Aligns effect/non-effect alleles to a common reference strand (flipping $\beta \to -\beta$ and $\text{EAF} \to 1 - \text{EAF}$ for swapped alleles, resolving unambiguous strand flips, dropping incompatible or EAF-ambiguous palindromic SNPs, and flagging cross-cohort EAF range $> \text{eaf\_diff\_max}$).
- **Fixed & Random Effects**:
  - `FE_IVW`: inverse-variance weighted $\hat{\beta}_{\text{FE}}$ with Cochran's $Q$ and $I^2 = \max\left(0, \frac{Q - (k - 1)}{Q}\right) \times 100\%$.
  - `RE_DL`: DerSimonian-Laird $\tau^2_{\text{DL}} = \max\left(0, \frac{Q - (k - 1)}{C}\right)$ random-effects estimator.

## 5. Polygenic Risk Score Evaluation (`gwas_prs_eval`)

- **Nested Model Comparison**: Compares **Clinical Covariates Only** vs **Covariates + `prs_z`** using $\Delta R^2$ (OLS), Nagelkerke pseudo-$\Delta R^2$ (logistic / Cox), Lee et al. (2012) liability $R^2$ when `pop_prev` is provided, likelihood-ratio test $P$, and $\Delta\text{AUC}$ or Harrell's $\Delta\text{C-index}$.
- **Quantile Risk Stratification**: Estimates covariate-adjusted OR / HR (or mean difference) with 95% CIs for `Bottom 10%`, `Top 20%`, `Top 10%`, and `Top 5%` relative to `Middle 40-60%`.
