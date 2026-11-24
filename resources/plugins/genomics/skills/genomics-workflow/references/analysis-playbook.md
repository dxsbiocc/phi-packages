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
