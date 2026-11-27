---
name: germline-workflow
description: >-
  Clinical rare-disease germline variant interpretation and pharmacogenomics
  workflow: trio/singleton Mendelian inheritance filtering (de novo, AR
  homozygous, phased in-trans vs unphased compound het, AD, X-linked); ACMG/AMP
  2015 + ClinGen Pejaver (2022) + Tavtigian (2020) Bayesian point-scale scoring
  with ACMG SF v3.2 flagging; IC-weighted HPO phenotype ranking; CPIC
  star-allele diplotypes, Activity Scores, and Level A dosing alerts; and
  ClinGen germline CNV/SV dosage scoring with trans SNV+DEL compound hit
  detection. Use for rare disease, germline VCF, trio, de novo, compound het,
  ACMG, PVS1, HPO, Exomiser, CPIC, PGx, star allele, CYP2C19, CYP2D6, DPYD,
  germline CNV, 罕见病, 遗传病, 家系分析, 致病性评级, 药物基因组, 星号等位基因.
phi:
  environment: plugin:germ
  attachTo:
    - ClinicalGermline
  scripts:
    - name: vcf_filter
      description: >-
        Filter single-sample or trio annotated VCF/VCF.gz (INFO/CSQ or INFO/ANN)
        or variant TSV/CSV by read depth, genotype quality, coding/splice
        consequence, and population AF, evaluating Mendelian inheritance fit
        (de novo, AR homozygous, phased in-trans vs unphased compound het, AD,
        and X-linked).
      run: [python, ./scripts/germ_vcf_filter.py]
      args:
        type: object
        additionalProperties: false
        required: [variants, output_dir]
        properties:
          variants:
            type: string
            format: input-path
            description: Multi-sample or single-sample annotated VCF (.vcf/.vcf.gz) or VEP/SnpEff TSV/CSV variant table.
          output_dir:
            type: string
            format: project-path
            description: Output directory for filtered germline candidates TSV, inheritance gene summary TSV, and filtering funnel figure.
          proband_id:
            type: string
            description: Proband sample ID (auto-picked if the input contains a single sample).
          mother_id:
            type: string
            description: Mother sample ID for trio de novo and trans compound-heterozygous phasing.
          father_id:
            type: string
            description: Father sample ID for trio de novo and trans compound-heterozygous phasing.
          proband_sex:
            type: string
            enum: [unknown, male, female]
            description: Proband biological sex for X-linked hemizygous vs heterozygous evaluation (default unknown).
          max_af_dominant:
            type: number
            description: Maximum population allele frequency for dominant/de novo candidates (default 0.0001).
          max_af_recessive:
            type: number
            description: Maximum population allele frequency for recessive/X-linked candidates (default 0.01).
          min_dp:
            type: integer
            description: Minimum genotype read depth DP (default 10).
          min_gq:
            type: integer
            description: Minimum genotype quality GQ (default 20).
      approval: write
      output: ./schemas/vcf-filter-result.json
      timeoutSeconds: 600
    - name: acmg_prioritize
      description: >-
        Classify rare germline variants under ACMG/AMP 2015 + ClinGen Pejaver
        (2022) + Tavtigian (2020) Bayesian point-scale rules (PVS1, PS1/PS5,
        PS2, PM2_Supporting, PM3, PM4, PP3, BA1, BS1, BP4, BP7) and flag
        ACMG SF v3.2 actionable secondary-finding genes.
      run: [python, ./scripts/germ_acmg_prioritize.py]
      args:
        type: object
        additionalProperties: false
        required: [variants, output_dir]
        properties:
          variants:
            type: string
            format: input-path
            description: Variant TSV/CSV from germ_vcf_filter or annotated variant table with gene, consequence, and gnomad_af.
          output_dir:
            type: string
            format: project-path
            description: Output directory for ACMG classified variants TSV and evidence summary figure.
          lof_genes:
            type: string
            format: input-path
            description: Optional TSV/TXT of known haploinsufficient / LoF-mechanism genes (augments built-in 80+ gene baseline).
      approval: write
      output: ./schemas/acmg-prioritize-result.json
      timeoutSeconds: 600
    - name: phenotype_rank
      description: >-
        Compute Information-Content-weighted HPO phenotypic similarity between
        patient HPO terms/symptoms and rare-disease genes, integrating ACMG
        points and inheritance fit into an Exomiser-inspired diagnostic score
        (0..100).
      run: [python, ./scripts/germ_phenotype_rank.py]
      args:
        type: object
        additionalProperties: false
        required: [variants, hpo_terms, output_dir]
        properties:
          variants:
            type: string
            format: input-path
            description: Variant TSV/CSV from germ_acmg_prioritize or germ_vcf_filter.
          hpo_terms:
            type: array
            items:
              type: string
            description: Patient HPO term IDs (e.g. HP:0001250, HP:0001263) or symptom keywords (e.g. seizures, cardiomyopathy).
          output_dir:
            type: string
            format: project-path
            description: Output directory for diagnostic gene ranking TSV, phenotype-matched variants TSV, and prioritization plot.
          hpo_annotations:
            type: string
            format: input-path
            description: Optional HPO genes_to_phenotype.txt TSV augmenting the built-in multi-system rare-disease HPO mapping.
      approval: write
      output: ./schemas/phenotype-rank-result.json
      timeoutSeconds: 600
    - name: pgx_star_allele
      description: >-
        Resolve CPIC star-allele diplotypes, Activity Scores, metabolizer and
        transporter phenotypes, and CPIC Level A clinical dosing alerts across
        CYP2C19, CYP2D6, CYP2C9, VKORC1, SLCO1B1, DPYD, TPMT, NUDT15, and
        UGT1A1 from pre-called diplotypes or rsID/VCF genotypes.
      run: [python, ./scripts/germ_pgx_star_allele.py]
      args:
        type: object
        additionalProperties: false
        required: [genotypes, output_dir]
        properties:
          genotypes:
            type: string
            format: input-path
            description: Germline VCF/VCF.gz or TSV/CSV with star-allele diplotypes (gene, diplotype) or rsID/position genotypes.
          output_dir:
            type: string
            format: project-path
            description: Output directory for CPIC diplotype phenotypes TSV, drug recommendations TSV, and actionability overview figure.
      approval: write
      output: ./schemas/pgx-star-allele-result.json
      timeoutSeconds: 600
    - name: cnv_sv_integrate
      description: >-
        Score germline CNV/SV calls (DEL, DUP, INV, BND) using ClinGen/ACMG
        2020 Riggs dosage-sensitivity criteria and cross-reference rare coding
        SNVs to uncover trans-acting SNV + exon-DEL compound heterozygotes.
      run: [python, ./scripts/germ_cnv_sv_integrate.py]
      args:
        type: object
        additionalProperties: false
        required: [sv_calls, output_dir]
        properties:
          sv_calls:
            type: string
            format: input-path
            description: Germline CNV/SV BED, TSV/CSV, or VCF from Manta, CNVkit, ExomeDepth, or GATK-gCNV.
          output_dir:
            type: string
            format: project-path
            description: Output directory for classified germline CNV/SV TSV, SNV+CNV compound hits TSV, and dosage summary figure.
          snv_candidates:
            type: string
            format: input-path
            description: Optional rare SNV/indel candidates TSV from germ_vcf_filter or germ_acmg_prioritize to detect SNV + DEL compound heterozygotes.
          gene_bed:
            type: string
            format: input-path
            description: Optional exon/gene coordinate BED/TSV augmenting built-in ClinGen dosage-sensitive gene intervals.
      approval: write
      output: ./schemas/cnv-sv-integrate-result.json
      timeoutSeconds: 600
---

# Clinical Germline & PGx Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for the Mendelian phasing rules, Tavtigian Bayesian point-scale equations, Pejaver calibrated computational cutoffs, IC-weighted HPO similarity formula, CPIC Level A guidelines, and ClinGen CNV/SV dosage scoring details.

## Tool Overview

```
VCF (CSQ/ANN) / VEP TSV ──germ_vcf_filter──────▶ Mendelian inheritance candidates + gene summary
Candidates TSV          ──germ_acmg_prioritize─▶ ACMG/AMP 2015 + Tavtigian points + SF v3.2 flags
ACMG/filtered TSV + HPO ──germ_phenotype_rank──▶ IC-weighted HPO sim + diagnostic score (0..100)
PGx VCF / diplotype TSV ──germ_pgx_star_allele─▶ CPIC diplotypes, Activity Scores, Level A dosing
CNV/SV BED/VCF (+ SNVs) ──germ_cnv_sv_integrate▶ Riggs dosage tiers + SNV+DEL compound hits
```

## Hard Red Lines

1. **Enforce trio phasing for `compound_heterozygous` & `PS2`**: only assign phased *in trans* `compound_heterozygous` (eligible for `PM3`) or `de_novo` (`PS2`) when parental genotypes verify the phase/origin; label singleton multi-hit genes `unphased_candidate_comp_het`.
2. **Use calibrated single-source computational evidence**: apply Pejaver (2022) REVEL/CADD/SpliceAI thresholds without double-counting correlated predictors, and report `acmg_points` and Bayesian `posterior_prob` ($P_0 = 0.10, O_{\text{PVS}} = 350$) transparently.
3. **Separate primary diagnostic ranking from ACMG SF v3.2 secondary findings**: use `germ_phenotype_rank` with patient HPO terms for diagnostic prioritization, while independently surfacing actionable `acmg_sf_reportable` findings.
4. **Always check SNV + CNV deletion compound hits in recessive phenotypes**: pass `--snv_candidates` into `germ_cnv_sv_integrate` whenever CNV/SV calls are available.
5. **Escalate CPIC `contraindicated_or_avoid` and `dose_reduction_required` alerts**: never bury `DPYD`, `CYP2C19`, `CYP2D6`, `TPMT`, or `NUDT15` high-risk phenotypes in unactionable summaries.
