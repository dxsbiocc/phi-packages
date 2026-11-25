# DNA Methylation Playbook

## 1. Assay-Specific Threshold Starting Points

| Modality | Recommended `min_cov` | Typical `delta_beta` | Notes |
| --- | --- | --- | --- |
| WGBS / EM-seq (tissue/cell line) | `10` (`5` for low-pass) | `0.15` – `0.20` | Merge into 500 bp DMRs (`min_cpgs=3`) |
| RRBS | `10` | `0.10` – `0.15` | Enriched in CpG islands & promoters |
| Illumina EPIC / 450K (tumor vs normal) | N/A (array $\beta$) | `0.15` – `0.20` | Supply `--probe_coords` for DMR calling |
| Peripheral blood EWAS | N/A (array $\beta$) | `0.03` – `0.05` | Always adjust `--covariates` for age, sex, and leukocyte composition |

## 2. Why $M$-values + $\Delta\beta$?

- $\beta = \frac{M}{M + U}$ is directly interpretable ($0.20 = 20\%$ methylation difference), while $M = \log_2\left(\frac{\beta}{1-\beta}\right)$ is homoscedastic for Gaussian linear models and empirical-Bayes variance moderation (Du et al., *BMC Bioinformatics* 2010).
- `meth_dmr` clips $\beta \in [10^{-4}, 1 - 10^{-4}]$ before computing $M$-values, fits Smyth's empirical-Bayes moderated $t$-test (with multi-covariate design matrices when `--covariates` are passed), and requires BOTH `fdr < fdr_cutoff` AND `|delta_beta| >= delta_beta` for confident DMCs.

## 3. Coordinate-Based DMR Aggregation

- Adjacent CpGs within `max_dist` (default `500` bp) on the same chromosome that share the same direction of change (`sign(delta_beta)`) and pass the seed significance threshold (`pvalue < 0.05`) are chained into candidate regions.
- Regions containing at least `min_cpgs` (default `3`) CpGs receive a regional combined significance via Stouffer's Z method and regional BH FDR, prioritizing spatially coherent epigenetic switches over noisy single-cytosine fluctuations.

## 4. Genomic Context Interpretation

- **Promoter CpG Island (`Island` / `Shore` overlapping `Promoter`)**: hypermethylation ($\Delta\beta > 0$) classically silences transcription; hypomethylation ($\Delta\beta < 0$) marks enhancer/promoter activation or oncogene derepression.
- **Gene Body (`GeneBody`)**: DNA methylation in transcribed exons/introns is deposited by DNMT3B via H3K36me3 and frequently correlates **positively** with expression and splicing fidelity.
- **Distal Intergenic / OpenSea**: global hypomethylation across `OpenSea` is a hallmark of cancer genomic instability and aging.
