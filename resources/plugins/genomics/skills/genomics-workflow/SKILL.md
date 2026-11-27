---
name: genomics-workflow
description: >-
  Downstream somatic cancer genomics workflow in Python: convert annotated
  somatic VCFs (VEP CSQ, SnpEff ANN, ANNOVAR) into GDC MAF tables; validate MAF
  tables, clinical metadata, and CNV segments; summarize cohort mutations,
  Ti/Tv spectra, and oncoplots; compute FoCR TMB and MATH heterogeneity; test
  co-occurrence/exclusivity and 10 Sanchez-Vega pathways; refit SBS96 signatures
  via NNLS; quantify CNV FGA; score scarHRD (LOH/TAI/LST) and MSI biomarker
  strata; and estimate CCF clonal vs subclonal architecture. Use for genomics,
  基因组, 肿瘤基因组, 体细胞突变, VCF转MAF, vcf2maf, MAF, maftools, oncoplot,
  TMB, 突变负荷, MATH, mutational signature, SBS96, COSMIC, CNV, FGA, HRD, LOH,
  TAI, LST, MSI, clonality, CCF, subclonal, 同源重组缺陷, 微卫星不稳定, 克隆演化.
phi:
  environment: plugin:genomics
  attachTo:
    - Genomics
  scripts:
    - name: vcf2maf
      description: >-
        Convert an annotated somatic VCF (.vcf or .vcf.gz from nf-core/sarek,
        VEP CSQ, SnpEff ANN, or ANNOVAR) into a GDC-compliant MAF TSV with
        severity-prioritized transcript selection, t_depth/t_ref_count/t_alt_count/t_vaf
        extraction, and optional 3-mer Context extraction from a reference FASTA
        for downstream SBS96 signature refitting.
      run: [python, ./scripts/gen_vcf2maf.py]
      args:
        type: object
        additionalProperties: false
        required: [vcf, output_maf]
        properties:
          vcf:
            type: string
            format: input-path
            description: Annotated somatic VCF file (.vcf or .vcf.gz).
          output_maf:
            type: string
            format: project-path
            description: Output GDC MAF TSV path (.maf or .tsv).
          tumor_sample:
            type: string
            description: Optional explicit tumor sample column or barcode.
          normal_sample:
            type: string
            description: Optional matched normal sample column name in the VCF.
          include_filtered:
            type: boolean
            description: Retain non-PASS VCF records instead of filtering to PASS / '.' only (default false).
          ref_fasta:
            type: string
            format: input-path
            description: Optional reference genome FASTA (.fa/.fasta) to populate 3-mer Context for SBS96.
          ncbi_build:
            type: string
            description: Genome build identifier written to NCBI_Build (default GRCh38).
      approval: write
      output: ./schemas/vcf2maf-result.json
      timeoutSeconds: 600
    - name: inspect
      description: >-
        Inspect a MAF (.maf/.maf.gz/.tsv) somatic mutation table, optional
        clinical metadata table, and optional CBS/CNVkit .seg file: validate
        mandatory GDC MAF fields, summarize non-synonymous vs silent counts,
        check VAF and depth availability, verify clinical sample overlap, and
        flag quality or formatting issues. Call first before downstream steps.
      run: [python, ./scripts/gen_inspect.py]
      args:
        type: object
        additionalProperties: false
        required: [maf]
        properties:
          maf:
            type: string
            format: input-path
            description: Somatic MAF file (.maf, .maf.gz, or .tsv).
          clinical:
            type: string
            format: input-path
            description: Optional clinical metadata TSV/CSV.
          sample_col:
            type: string
            description: Clinical sample ID column (default Tumor_Sample_Barcode or first column).
          vaf_col:
            type: string
            description: Optional explicit VAF column name in the MAF.
          seg:
            type: string
            format: input-path
            description: Optional CBS/CNVkit copy-number segment file (.seg/.tsv).
      approval: read
      output: ./schemas/inspect-result.json
      timeoutSeconds: 300
    - name: summary
      description: >-
        Summarize cohort somatic mutations (non-synonymous variant
        classifications, per-sample burden, top mutated genes), compute the
        6-class SNV substitution spectrum and Ti/Tv ratio, and render a MAF
        summary dashboard and waterfall oncoplot (tracking Multi_Hit genes).
      run: [python, ./scripts/gen_summary.py]
      args:
        type: object
        additionalProperties: false
        required: [maf, output_dir]
        properties:
          maf:
            type: string
            format: input-path
            description: Somatic MAF file (.maf, .maf.gz, or .tsv).
          top_genes:
            type: integer
            minimum: 5
            description: Number of top mutated genes to show in the oncoplot (default 20).
          include_silent:
            type: boolean
            description: Include Silent/non-coding variants in Ti/Tv spectrum (default true).
          output_dir:
            type: string
            format: project-path
            description: Directory for gene/sample/TiTv tables and summary/oncoplot PNGs.
      approval: write
      output: ./schemas/summary-result.json
      timeoutSeconds: 600
    - name: tmb_ith
      description: >-
        Compute FoCR-compliant Tumor Mutational Burden (TMB = coding
        non-synonymous mutations / capture_size_mb after VAF filtering) and
        Mroz/Rocco Mutant-Allele Tumor Heterogeneity (MATH = 100 * MAD(VAF) /
        median(VAF)) per sample, optionally excluding variants inside
        CNV-altered segments (.seg) and testing group differences across a
        clinical variable.
      run: [python, ./scripts/gen_tmb_ith.py]
      args:
        type: object
        additionalProperties: false
        required: [maf, capture_size_mb, output_dir]
        properties:
          maf:
            type: string
            format: input-path
            description: Somatic MAF file (.maf, .maf.gz, or .tsv).
          capture_size_mb:
            type: number
            description: Explicit callable coding territory in Mb (e.g. 38.0 for WES, 1.2 for a targeted panel).
          min_vaf:
            type: number
            description: Minimum VAF filter for TMB and MATH (default 0.05; ignored if VAF absent).
          vaf_col:
            type: string
            description: Optional explicit VAF column name in the MAF.
          seg:
            type: string
            format: input-path
            description: Optional CBS/CNVkit .seg file to exclude CNV-altered loci from MATH.
          cnv_seg_cutoff:
            type: number
            description: Absolute Segment_Mean log2 ratio threshold defining CNV-altered regions (default 0.2).
          clinical:
            type: string
            format: input-path
            description: Optional clinical metadata TSV/CSV for group comparison.
          sample_col:
            type: string
            description: Clinical sample ID column (default Tumor_Sample_Barcode or first column).
          group_col:
            type: string
            description: Optional clinical column for Mann-Whitney / Kruskal-Wallis comparison of TMB and MATH.
          output_dir:
            type: string
            format: project-path
            description: Directory for per-sample TMB/MATH table and distribution figures.
      approval: write
      output: ./schemas/tmb-ith-result.json
      timeoutSeconds: 600
    - name: interactions
      description: >-
        Run pairwise Fisher's exact tests for somatic co-occurrence and mutual
        exclusivity (with BH FDR), summarize alterations across the 10
        canonical Sanchez-Vega (Cell 2018) oncogenic signaling pathways, and
        optionally test differential mutation frequencies between two clinical
        groups or two MAF cohorts (Fisher + odds ratio + forest plot).
      run: [python, ./scripts/gen_interactions.py]
      args:
        type: object
        additionalProperties: false
        required: [maf, output_dir]
        properties:
          maf:
            type: string
            format: input-path
            description: Primary somatic MAF file (.maf, .maf.gz, or .tsv).
          compare_maf:
            type: string
            format: input-path
            description: Optional second cohort MAF for two-cohort differential mutation comparison.
          clinical:
            type: string
            format: input-path
            description: Optional clinical metadata TSV/CSV for within-cohort group comparison.
          sample_col:
            type: string
            description: Clinical sample ID column.
          group_col:
            type: string
            description: Clinical grouping column for differential mutation testing.
          case:
            type: string
            description: Case group value in group_col (required when comparing clinical groups).
          control:
            type: string
            description: Control group value in group_col (required when comparing clinical groups).
          top_genes:
            type: integer
            minimum: 5
            description: Top mutated genes for pairwise interaction matrix (default 20).
          min_mut:
            type: integer
            minimum: 2
            description: Minimum mutated samples required to test a gene in cohort comparison (default 3).
          fdr:
            type: number
            description: BH FDR cutoff for highlighting significant pairs/genes (default 0.05).
          output_dir:
            type: string
            format: project-path
            description: Directory for interactions, pathways, cohort comparison tables, and plots.
      approval: write
      output: ./schemas/interactions-result.json
      timeoutSeconds: 600
    - name: signatures
      description: >-
        Extract the 96-channel SBS trinucleotide mutational catalog from a MAF
        carrying a Context column (or user-supplied reference FASTA) and refit
        per-sample exposures against a reference signature matrix (e.g. COSMIC
        v3.4 SBS TSV) using scipy.optimize.nnls with cosine similarity
        reconstruction diagnostics. When trinucleotide context is unavailable,
        reports the 6-class substitution spectrum and states the gap explicitly.
      run: [python, ./scripts/gen_signatures.py]
      args:
        type: object
        additionalProperties: false
        required: [maf, output_dir]
        properties:
          maf:
            type: string
            format: input-path
            description: Somatic MAF file (.maf, .maf.gz, or .tsv).
          context_col:
            type: string
            description: Optional MAF column containing trinucleotide context (e.g. Context, Triplet).
          ref_signatures:
            type: string
            format: input-path
            description: Optional reference SBS96 signature matrix TSV/CSV (Type + signature columns, e.g. COSMIC v3.4).
          min_snvs:
            type: integer
            minimum: 1
            description: Minimum SNVs per sample required for NNLS signature refitting (default 20).
          output_dir:
            type: string
            format: project-path
            description: Directory for SBS96 catalog, exposure matrix, cosine similarity table, and plots.
      approval: write
      output: ./schemas/signatures-result.json
      timeoutSeconds: 600
    - name: cnv
      description: >-
        Analyze CBS or CNVkit copy-number segment files (.seg/.tsv): compute
        per-sample Fraction Genome Altered (FGA), chromosome-level gain/loss
        burden, recurrent altered chromosomes/segments, and optional clinical
        group FGA comparisons.
      run: [python, ./scripts/gen_cnv.py]
      args:
        type: object
        additionalProperties: false
        required: [seg, output_dir]
        properties:
          seg:
            type: string
            format: input-path
            description: Copy-number segment file (.seg/.tsv with Sample, Chromosome, Start, End, Segment_Mean).
          gain_threshold:
            type: number
            description: Log2 ratio threshold for copy-number gain (default 0.2).
          loss_threshold:
            type: number
            description: Log2 ratio threshold for copy-number loss (default -0.2).
          amp_threshold:
            type: number
            description: Log2 ratio threshold for high-level amplification (default 0.8).
          del_threshold:
            type: number
            description: Log2 ratio threshold for deep deletion (default -0.8).
          clinical:
            type: string
            format: input-path
            description: Optional clinical metadata TSV/CSV for FGA group comparison.
          sample_col:
            type: string
            description: Clinical sample ID column.
          group_col:
            type: string
            description: Clinical grouping column for FGA comparison.
          output_dir:
            type: string
            format: project-path
            description: Directory for per-sample FGA table, chromosome summary, and CNV figures.
      approval: write
      output: ./schemas/cnv-result.json
      timeoutSeconds: 600
    - name: hrd_msi
      description: >-
        Compute scarHRD/Marquard genomic scar Homologous Recombination Deficiency
        scores (HRD-LOH, TAI, LST, and hrd_sum on hg38/hg19 autosomes) from
        allele-specific CNV segments (ASCAT/Sequenza/FACETS/PURPLE) and integrate
        MSIsensor-pro/Mantis MSI tables or somatic MAF indel/frameshift surrogates
        to classify combined HRD/MSI biomarker strata.
      run: [python, ./scripts/gen_hrd_msi.py]
      args:
        type: object
        additionalProperties: false
        required: [output_dir]
        properties:
          segments:
            type: string
            format: input-path
            description: Optional allele-specific CNV segment TSV/CSV (sample, chrom, start, end, nMajor/nMinor or total_cn/minor_cn).
          msi_table:
            type: string
            format: input-path
            description: Optional MSIsensor-pro or Mantis MSI summary TSV/CSV.
          maf:
            type: string
            format: input-path
            description: Optional somatic MAF/TSV for frameshift indel rate and indel-to-SNV MSI surrogate.
          genome:
            type: string
            enum: [hg38, hg19]
            description: Reference genome assembly for chromosome lengths and centromeres (default hg38).
          hrd_cutoff:
            type: integer
            description: HRD sum threshold for HRD-positive status (default 42).
          msi_pct_cutoff:
            type: number
            description: Unstable microsatellite percentage threshold for MSI-H (default 3.5).
          output_dir:
            type: string
            format: project-path
            description: Directory for hrd_msi_summary.tsv and hrd_msi_landscape.png.
      approval: write
      output: ./schemas/hrd-msi-result.json
      timeoutSeconds: 600
    - name: clonality
      description: >-
        Estimate McGranahan/PyClone-inspired Cancer Cell Fraction (CCF) with
        Binomial multiplicity MLE and exact Clopper-Pearson 95% confidence
        intervals from somatic variant counts, sample purity, and local copy
        number; classify clonal vs subclonal driver mutations; and test
        Williams et al. (2016) neutral tail linearity M(f) ~ 1/f.
      run: [python, ./scripts/gen_clonality.py]
      args:
        type: object
        additionalProperties: false
        required: [maf, output_dir]
        properties:
          maf:
            type: string
            format: input-path
            description: Somatic MAF/TSV with Hugo_Symbol, Tumor_Sample_Barcode, Chromosome, Start_Position, and counts or VAF+depth.
          output_dir:
            type: string
            format: project-path
            description: Directory for clonality_variants.tsv, clonality_sample_summary.tsv, and ccf_clonal_architecture.png.
          purity_table:
            type: string
            format: input-path
            description: Optional sample purity/ploidy TSV/CSV (sample, purity, ploidy).
          segments:
            type: string
            format: input-path
            description: Optional CNV segment TSV/CSV (sample, chrom, start, end, total_cn or nMajor+nMinor).
          default_purity:
            type: number
            description: Fallback tumor purity when purity_table is omitted or sample is absent (default 0.70).
          default_cn:
            type: integer
            description: Fallback local total copy number when segments are omitted or locus uncovered (default 2).
          clonal_ccf_cutoff:
            type: number
            description: CCF threshold for clonal classification alongside 95% CI upper bound >= 0.95 (default 0.85).
      approval: write
      output: ./schemas/clonality-result.json
      timeoutSeconds: 600
---

# Genomics — somatic cancer genomics downstream chain

Read [references/analysis-playbook.md](references/analysis-playbook.md) before
interpreting TMB, MATH heterogeneity, oncogenic pathways, SBS96 signatures,
HRD/MSI stratification, or CCF clonality.

## Chain

```
FASTQ ──(nf-core/sarek / Mutect2 / VEP / vcf2maf)──▶ MAF + CNV .seg / ASCAT + MSI
MAF (+ clinical + seg) ──gen_inspect───────▶ validate GDC MAF fields, VAF, sample overlap
                       ──gen_summary───────▶ cohort summary, Ti/Tv spectrum & oncoplot
                       ──gen_tmb_ith───────▶ FoCR TMB (explicit Mb) & MATH heterogeneity
                       ──gen_interactions──▶ co-occurrence/exclusivity, 10 pathways & cohort diff
                       ──gen_signatures────▶ SBS96 catalog & COSMIC NNLS refitting
                       ──gen_clonality─────▶ CCF multiplicity MLE, clonal/subclonal drivers & 1/f tail
CNV .seg (+ clinical)  ──gen_cnv───────────▶ per-sample FGA & recurrent gain/loss summary
ASCAT/MSI/MAF          ──gen_hrd_msi───────▶ scarHRD (LOH/TAI/LST), MSI & biomarker strata
```

## Hard red lines

1. **Mandatory GDC MAF fields**: verify `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Classification`, `Variant_Type`, and `Tumor_Sample_Barcode` via `gen_inspect`.
2. **Never invent `capture_size_mb` for TMB**: WES (~35–50 Mb) and targeted panels (~0.5–3 Mb) yield completely different denominators; always pass the assay's actual callable coding territory.
3. **Separate oncogenicity from clinical actionability**: pathway hits are biological driver evidence, not Tier I clinical therapy recommendations.
4. **Never fabricate SBS96 exposures**: without trinucleotide context and a reference signature matrix, report the 6-class Ti/Tv fallback and the exact gap.
5. **Control for CNV in MATH**: pass `--seg` to `gen_tmb_ith` whenever copy-number segments are available.
6. **Disclose HRD/MSI modality limits**: allele-specific segments are required for scarHRD (`LOH`/`TAI`/`LST`), and MAF indel ratios without an MSIsensor-pro/Mantis table are proxy surrogates only.
7. **State purity/CN inputs in CCF clonality**: report whether purity and local copy number came from measured tables or `default_purity`/`default_cn` fallbacks.
