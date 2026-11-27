# Changelog — statistical-genetics

## 0.1.0 (2026-11-27)

Initial release of the `statistical-genetics` plugin (`toolPrefix: gwas`, agent `StatisticalGenetics`, private locked environment `plugin:gwas`).

- `gwas_sumstats_harmonize`: Multi-caller GWAS summary-statistics harmonization (PLINK2, REGENIE, SAIGE, BOLT-LMM, GWAS Catalog), `OR -> ln(OR)` & `LOG10P -> P` conversion, IUPAC & `SE > 0` QC, $\lambda_{\text{GC}}$ / $\lambda_{1000}$ inflation metrics, physical-window lead locus extraction, and Manhattan/QQ plotting.
- `gwas_ldsc_interpret`: Heteroscedasticity-weighted univariate and bivariate LD Score Regression with genomic block jackknife SEs, Lee et al. (2011) liability-scale $h^2$, LDSC intercept & attenuation ratio, and cross-trait genetic correlation $r_g$.
- `gwas_finemap`: SuSiE-RSS Iterative Bayesian Stepwise Selection with regularized LD $\tilde{R} = (1-\lambda)R + \lambda I$, purity-filtered 95% credible sets, and Wakefield ABF $L=1$ fallback when no LD matrix is supplied.
- `gwas_meta_qc`: Cross-cohort allele alignment, EAF outlier flagging, Fixed-Effect IVW and DerSimonian-Laird Random-Effects meta-analysis with Cochran's $Q$, $I^2$, and direction strings.
- `gwas_prs_eval`: Standardized PRS evaluation comparing Covariates Only vs Covariates + PRS across linear, logistic, and Cox PH survival models with $\Delta R^2$, liability $R^2$, $\Delta\text{AUC}$ / $\Delta\text{C-index}$, and quantile risk stratification.
