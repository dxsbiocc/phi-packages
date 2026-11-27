# Somatic Genomics Downstream Analysis Playbook

## 1. MAF Validation (`gen_inspect`)
- Per the NCI GDC MAF specification and `maftools` (Mayakonda et al., Genome Research 2018), nine fields are mandatory: `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Classification`, `Variant_Type`, and `Tumor_Sample_Barcode`.
- Standard non-synonymous variant classifications (`vc_nonSyn`): `Frame_Shift_Del`, `Frame_Shift_Ins`, `In_Frame_Del`, `In_Frame_Ins`, `Missense_Mutation`, `Nonsense_Mutation`, `Nonstop_Mutation`, `Splice_Site`, and `Translation_Start_Site`.
- VAF may appear directly (`t_vaf`, `VAF`, `i_TumorVAF_WU`) or be computed from `t_alt_count / t_depth` (or `t_alt_count / (t_alt_count + t_ref_count)`). Percentage-scale VAFs (> 1) are automatically normalized to [0, 1].

## 2. Cohort Landscape (`gen_summary`)
- **Ti/Tv**: transitions (`C>T`, `T>C`) vs transversions (`C>A`, `C>G`, `T>A`, `T>G`), folded onto pyrimidine reference bases (`C`, `T`) using reverse complements for `G`/`A` references.
- **Oncoplot**: ranks top mutated genes by fraction of mutated samples; multiple non-synonymous hits in the same gene and sample are labeled `Multi_Hit`.

## 3. TMB & Intratumor Heterogeneity (`gen_tmb_ith`)
- **TMB (FoCR Harmonization Project, JITC 2020)**: `TMB = non_synonymous_coding_mutations / capture_size_mb`, filtering low-VAF artifacts (`min_vaf = 0.05` default when VAF is present). Never compare raw WES and targeted-panel TMB without noting panel size and filtering parity.
- **MATH (Mroz & Rocco)**: `MATH = 100 * MAD(VAF) / median(VAF)`, where `MAD = median(|VAF_i - median(VAF)|)`. Higher MATH reflects broader subclonal VAF dispersion. Pass `--seg` to exclude variants in CNV-altered segments (`|Segment_Mean| >= 0.2`) so copy-number shifts do not inflate MATH.

## 4. Somatic Interactions & Oncogenic Pathways (`gen_interactions`)
- **Co-occurrence / Mutual Exclusivity**: pairwise 2×2 Fisher's exact test across the top mutated genes with Benjamini-Hochberg FDR; odds ratio > 1 indicates co-occurrence, < 1 indicates mutual exclusivity.
- **10 Canonical Oncogenic Signaling Pathways (Sanchez-Vega et al., Cell 2018)**: Cell Cycle, Hippo, Myc, Notch, Nrf2, PI3K/Akt, RTK-RAS, TGF-β, p53, and Wnt/β-catenin.
- **Cohort Comparison**: gene-wise Fisher's exact test comparing mutation prevalence between two MAF cohorts or two clinical groups, with odds ratio, 95% confidence interval, and BH FDR.

## 5. Mutational Signatures (`gen_signatures`)
- **SBS96**: classifies each SNV into 96 pyrimidine-centered trinucleotide contexts `5'[Ref>Alt]3'`.
- **NNLS Refitting**: solves `min_{x >= 0} || W x - y ||_2` via `scipy.optimize.nnls` against a user-supplied reference signature matrix (e.g. COSMIC v3.4 SBS96), reporting both relative exposures and cosine similarity between observed and reconstructed 96-channel spectra (cosine < 0.8 indicates poor reconstruction or low mutation count).

## 6. Copy-Number Segments (`gen_cnv`)
- **Fraction Genome Altered (FGA)**: sum of segment lengths with `|Segment_Mean| >= gain/loss_threshold` divided by total profiled genomic length per sample.
- Distinguishes broad low-level gain/loss (`±0.2`) from high-level amplification/deep deletion (`±0.8`).

## 7. Genomic Scar HRD & MSI Stratification (`gen_hrd_msi`)
- **scarHRD / Marquard Three-Component Score (autosomes `1–22`, `hg38`/`hg19`)**:
  - **HRD-LOH** (Abkevich et al., 2012): contiguous LOH regions (`minor_cn == 0`, `total_cn >= 1`) strictly `> 15 Mb` that do not span an entire chromosome.
  - **TAI** (Birkbak et al., 2012): allelic-imbalance segments (`major_cn != minor_cn`) extending to within `2 Mb` of a telomere with length `>= 11 Mb`, excluding centromere-crossing and whole-arm/whole-chromosome events.
  - **LST** (Popova et al., 2012): chromosomal breaks between adjacent segments each `>= 10 Mb` (after smoothing segments `< 3 Mb`) separated by `<= 3 Mb` on the same chromosome arm (`p` or `q`).
  - **`hrd_sum = LOH + TAI + LST`**: clinical positivity evaluated at `hrd_cutoff = 42` (Myriad/Marquard threshold).
- **MSI & Combined Biomarker Strata**: integrates MSIsensor-pro/Mantis `%` unstable sites (`msi_pct >= 3.5` → `MSI-H`) or, when absent, somatic MAF `frameshift_indel_rate` and `indel_to_snv_ratio` (`MSI-H_proxy` surrogate) to stratify `HRD_pos_MSS` (PARP-inhibitor candidate), `MSI_H` (immunotherapy candidate), and `HR_proficient_MSS`.

## 8. Cancer Cell Fraction (CCF) & Clonal Architecture (`gen_clonality`)
- **Expected Clonal VAF & Multiplicity MLE (McGranahan et al. / PyClone model)**: for tumor purity $\rho$ and local total copy number $C_{\text{tot}}$, $f_{\text{exp}}(m) = \frac{\rho m}{\rho C_{\text{tot}} + 2(1-\rho)}$ for $m \in \{1, \dots, \max(1, C_{\text{tot}})\}$. Selects $\hat{m} = \arg\max_{m \ge 1} \text{Binom}(k_{\text{alt}} \mid N_{\text{depth}}, \min(0.999, f_{\text{exp}}(m)))$ when $k_{\text{alt}}/N \ge f_{\text{exp}}(1)$ (else $\hat{m} = 1$), and estimates $\widehat{\text{CCF}} = \min\left(1.5, \frac{\text{VAF} [\rho C_{\text{tot}} + 2(1-\rho)]}{\rho \hat{m}}\right)$ with exact Clopper-Pearson 95% confidence bounds.
- **Clonality Call**: **`clonal`** if `ccf_ci_high >= 0.95` or `ccf >= clonal_ccf_cutoff` (`0.85`), else **`subclonal`**.
- **Neutral Tail Test (Williams et al., Nature Genetics 2016)**: evaluates linearity of cumulative mutation count $M(f) \propto 1/f - 1/f_{\max}$ across $0.12 \le \text{VAF} \le 0.24$, reporting $R^2_{\text{neutral}}$.
