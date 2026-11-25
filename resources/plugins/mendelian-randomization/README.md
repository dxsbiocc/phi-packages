# mendelian-randomization · 0.1.0

Genetic epidemiology, GWAS summary-statistics QC, Two-Sample Mendelian
Randomization (MR), and Bayesian colocalization plugin (`toolPrefix: mr`, agent
`MendelianRandomization`, private locked Python 3.13 environment `plugin:mr`).

## Tools (`mr_*`)

1. `mr_gwas_qc` — standardize GWAS summary statistics (`BETA` or `OR -> log(OR)`,
   allele uppercase/indel cleaning), compute genomic inflation $\lambda_{\text{GC}}$,
   extract distance-clumped lead variants, and render Manhattan + QQ plots.
2. `mr_harmonize` — align exposure genetic instruments (`F >= 10`) against
   outcome GWAS alleles, resolve strand flips and palindromic `A/T`/`C/G` SNPs
   via EAF margins, and optionally apply Steiger variance-explained filtering.
3. `mr_estimate` — run multiplicative random-effects IVW, fixed-effects IVW (or
   Wald ratio for 1-SNP cis-pQTL/eQTL), MR-Egger regression + intercept
   pleiotropy test, Weighted Median, Simple/Weighted Mode, Cochran's $Q$ & $I^2$,
   residual outlier flags, Leave-One-Out sensitivity, and 4-panel MR plots.
4. `mr_coloc` — compute exact Wakefield (2009) / Giambartolomei (2014) Bayesian
   colocalization (`PP.H0`–`PP.H4`) across a shared cis locus with 95% credible
   set identification and regional mirror plots.
