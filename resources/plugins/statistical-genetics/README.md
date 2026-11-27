# statistical-genetics · 0.1.0

Statistical genetics, GWAS summary-statistics harmonization, LD Score Regression,
multi-causal SuSiE-RSS fine-mapping, cross-cohort meta-analysis, and polygenic
risk score (PRS) evaluation plugin (`toolPrefix: gwas`, agent
`StatisticalGenetics`) running in a locked MIT/BSD `conda-forge` environment
(`plugin:gwas`: Python 3.13 + `pandas`/`numpy`/`scipy`/`scikit-learn`/`matplotlib-base`).

All code in `plugin:gwas` is 100% original MIT-licensed Python.

## Division of Work: Upstream Wrappers, `gwas_*` Tools, and Downstream `mendelian-randomization` (`mr_*`)

| Stage / Task | Upstream `nf-core` Wrapper | Deterministic Local Tool (`plugin:gwas`) | Downstream `mendelian-randomization` (`plugin:mr`) |
| --- | --- | --- | --- |
| Genotype QC, PCA, MAF/HWE filtering & LD matrix computation | `nf-core/modules/plink2` | — | — |
| Whole-genome mixed-model association (quantitative & binary) | `nf-core/modules/regenie` | — | — |
| GREML SNP heritability, COJO conditional & fastGWA | `nf-core/modules/gcta` | — | — |
| Genotype phasing & imputation (large biobanks / reference panels) | `nf-core/modules/shapeit5`, `nf-core/modules/beagle5` | — | — |
| GWAS summary-stat harmonization, $\lambda_{\text{GC}}/\lambda_{1000}$, lead clumping | — | `gwas_sumstats_harmonize` | — |
| Univariate & cross-trait bivariate LD Score Regression ($h^2$, intercept, $r_g$) | — | `gwas_ldsc_interpret` | — |
| Multi-causal SuSiE-RSS / Wakefield ABF fine-mapping & 95% credible sets | — | `gwas_finemap` | — |
| Multi-cohort allele alignment, EAF QC, FE-IVW & RE-DL meta-analysis ($Q$, $I^2$) | — | `gwas_meta_qc` | — |
| PRS incremental $\Delta R^2$/AUC/C-index & quantile risk stratification | — | `gwas_prs_eval` | — |
| Two-sample Mendelian Randomization (IVW/Egger/Median) & cis-coloc (`PP.H4`) | — | — | `mr_harmonize`, `mr_estimate`, `mr_coloc` |

## Tools (`gwas_*`)

- `gwas_sumstats_harmonize`: Standardizes PLINK2, REGENIE, SAIGE, BOLT-LMM, and GWAS Catalog summary tables (`BETA`, `OR -> ln(OR)`, `LOG10P -> P`, `Z`, `N_EFF`), enforces IUPAC alleles & `SE > 0`, applies `min_maf`/`min_info` QC, computes $\lambda_{\text{GC}}$ and $\lambda_{1000}$, extracts distance-clumped lead loci, and renders a 2-panel Manhattan + QQ figure.
- `gwas_ldsc_interpret`: Runs Bulik-Sullivan (2015) heteroscedasticity-weighted univariate and bivariate LD Score Regression with genomic block jackknife SEs, estimating observed/liability $h^2$, LDSC intercept, attenuation ratio, and cross-trait genetic correlation $r_g$.
- `gwas_finemap`: Runs SuSiE-RSS Iterative Bayesian Stepwise Selection on summary $Z$-scores and regularized LD $\tilde{R} = (1-\lambda)R + \lambda I$ (or Wakefield ABF $L=1$ when `ld_matrix` is omitted), reporting per-SNP PIPs and purity-filtered 95% credible sets.
- `gwas_meta_qc`: Aligns alleles across $\ge 2$ GWAS cohorts (with strand-flip and EAF-discordance QC), runs Fixed-Effect IVW and DerSimonian-Laird Random-Effects meta-analysis, computes Cochran's $Q$, $I^2$, and direction strings, and renders heterogeneity/EAF diagnostic plots.
- `gwas_prs_eval`: Evaluates standardized `prs_z` across binary, quantitative, or time-to-event Cox PH outcomes, comparing Clinical Covariates Only vs Covariates + PRS ($\Delta R^2$, Nagelkerke/liability $\Delta R^2$, LRT $P$, $\Delta\text{AUC}$ / $\Delta\text{C-index}$) and estimating OR/HR across PRS quantiles vs the `Middle 40-60%` reference.
