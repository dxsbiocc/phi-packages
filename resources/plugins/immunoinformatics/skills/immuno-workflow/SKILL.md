---
name: immuno-workflow
description: >-
  Downstream TCR/BCR immune repertoire, HLA genotyping consensus, and neoantigen
  prioritization workflow: AIRR/10x/TRUST4/MiXCR repertoire harmonization and
  spectratype QC; Shannon/Pielou/Simpson/D50 diversity, clonal homeostasis space
  and V-J pairing; Morisita-Horn/Jaccard overlap, 4-mer motifs and Fisher
  longitudinal clonotype tracking; 2-field HLA consensus across OptiType,
  arcasHLA, seq2HLA and HLA-LA with Class I homozygosity and approximate HED;
  and neoantigen MHC IC50, agretopicity (DAI), anchor/foreignness and TPM/VAF
  prioritization. Use for TCR, BCR, repertoire, AIRR, TRUST4, MiXCR, HLA, HED,
  neoantigen, agretopicity, immunogenicity, 免疫组库, 克隆型, 新抗原, 人类白细胞抗原.
phi:
  environment: plugin:imm
  attachTo:
    - Immunoinformatics
  scripts:
    - name: repertoire_inspect
      description: >-
        Harmonize AIRR, 10x Cell Ranger V(D)J, TRUST4, and MiXCR repertoire
        tables into standardized AIRR-compatible TSV with per-sample QC metrics
        (productive %, out-of-frame/stop %, chain breakdown, CDR3 AA length &
        Kyte-Doolittle hydropathy) and spectratype QC figure.
      run: [python, ./scripts/imm_repertoire_inspect.py]
      args:
        type: object
        additionalProperties: false
        required: [inputs, output_dir]
        properties:
          inputs:
            type: array
            items:
              type: string
              format: input-path
            description: Input repertoire TSV/CSV files (AIRR, 10x filtered_contig_annotations.csv, TRUST4 report, or MiXCR clones).
          output_dir:
            type: string
            format: project-path
            description: Output directory for harmonized repertoire TSV, per-sample QC summary, and spectratype figure.
          format:
            type: string
            enum: [auto, airr, cellranger_vdj, trust4, mixcr]
            description: Input repertoire table format (default auto detects per file).
          chain_filter:
            type: string
            enum: [all, TCR, BCR, TRB, TRA, IGH]
            description: Receptor chain filter applied to harmonized repertoire output (default all).
          productive_only:
            type: boolean
            description: Retain only productive clonotypes in harmonized_repertoire.tsv after computing sample QC metrics (default false).
      approval: write
      output: ./schemas/repertoire-inspect-result.json
      timeoutSeconds: 600
    - name: diversity
      description: >-
        Compute per-sample repertoire diversity (Shannon entropy in bits and
        nats, Pielou clonality, Gini-Simpson, Inverse Simpson, D50), clonal
        homeostasis space tiers, V-J gene pairing usage matrix, and optional
        Mann-Whitney / Kruskal-Wallis + BH FDR group comparisons.
      run: [python, ./scripts/imm_diversity.py]
      args:
        type: object
        additionalProperties: false
        required: [repertoire, output_dir]
        properties:
          repertoire:
            type: string
            format: input-path
            description: Harmonized or raw repertoire TSV/CSV file.
          output_dir:
            type: string
            format: project-path
            description: Output directory for diversity metrics, V-J pairing matrix, optional group stats, and 2-panel homeostasis figure.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV containing sample_id/sample and grouping columns.
          group_col:
            type: string
            description: Metadata column name defining sample groups for statistical comparison.
          locus:
            type: string
            description: Target receptor locus to analyze (default TRB; pass ALL to pool loci).
      approval: write
      output: ./schemas/diversity-result.json
      timeoutSeconds: 600
    - name: clonotype_tracking
      description: >-
        Compute pairwise Morisita-Horn and Jaccard CDR3AA+V overlap matrices,
        public/shared clonotypes and enriched 4-mer CDR3 motifs, plus two-sided
        Fisher's exact test + BH FDR calling significantly expanded vs
        contracted clonotypes between baseline and follow-up samples.
      run: [python, ./scripts/imm_clonotype_tracking.py]
      args:
        type: object
        additionalProperties: false
        required: [repertoire, output_dir]
        properties:
          repertoire:
            type: string
            format: input-path
            description: Multi-sample harmonized or raw repertoire TSV/CSV file.
          output_dir:
            type: string
            format: project-path
            description: Output directory for overlap matrix, differential clonotype tracking table, and overlap figure.
          baseline_sample:
            type: string
            description: Baseline sample_id for longitudinal expansion/contraction testing (auto-picked if exactly 2 samples exist).
          followup_sample:
            type: string
            description: Follow-up sample_id for longitudinal expansion/contraction testing (auto-picked if exactly 2 samples exist).
          locus:
            type: string
            description: Target receptor locus to track (default TRB; pass ALL to pool loci).
          top_n:
            type: integer
            description: Number of top clonotypes highlighted in plots and summary (default 20).
      approval: write
      output: ./schemas/clonotype-tracking-result.json
      timeoutSeconds: 600
    - name: hla_summary
      description: >-
        Standardize OptiType, arcasHLA, seq2HLA, HLA-LA, and cohort HLA tables
        to 2-field (4-digit) alleles across Class I (HLA-A/B/C) and Class II
        (HLA-DRB1/DQB1/DPB1), flag Class I homozygosity, and compute approximate
        pocket-Grantham HLA Evolutionary Divergence (HED).
      run: [python, ./scripts/imm_hla_summary.py]
      args:
        type: object
        additionalProperties: false
        required: [hla_inputs, output_dir]
        properties:
          hla_inputs:
            type: array
            items:
              type: string
              format: input-path
            description: Input HLA typing files (OptiType TSV, arcasHLA JSON, seq2HLA genotype TSV, HLA-LA bestguess TXT, or clinical HLA TSV/CSV).
          output_dir:
            type: string
            format: project-path
            description: Output directory for consensus HLA genotypes TSV and zygosity/HED figure.
      approval: write
      output: ./schemas/hla-summary-result.json
      timeoutSeconds: 600
    - name: neoantigen_prioritize
      description: >-
        Prioritize candidate neoantigens from upstream MHC binding predictions
        or somatic peptide tables (with built-in anchor-motif PSSM IC50
        estimation), computing Differential Agretopicity Index (DAI), anchor vs
        TCR-contact mutation status, foreignness/hydropathy, TPM/VAF/CCF gates,
        and composite immunogenicity score (0..100).
      run: [python, ./scripts/imm_neoantigen_prioritize.py]
      args:
        type: object
        additionalProperties: false
        required: [candidates, output_dir]
        properties:
          candidates:
            type: string
            format: input-path
            description: Candidate TSV/CSV with upstream MHC predictions (mt_epitope, wt_epitope, hla_allele, ic50_mt, ic50_wt) or somatic peptides (gene, wt_peptide, mt_peptide, optional tpm/vaf/ccf).
          output_dir:
            type: string
            format: project-path
            description: Output directory for prioritized neoantigens TSV and neoantigen landscape figure.
          hla_alleles:
            type: array
            items:
              type: string
            description: Optional patient HLA alleles (e.g. HLA-A*02:01, HLA-B*07:02) for built-in anchor-motif PSSM IC50 estimation when ic50_mt is absent.
          ic50_cutoff_nm:
            type: number
            description: Binding affinity gate in nM for weak/strong binders (default 500.0).
          min_tpm:
            type: number
            description: Minimum transcript expression (TPM) gate when TPM is supplied (default 1.0).
      approval: write
      output: ./schemas/neoantigen-prioritize-result.json
      timeoutSeconds: 600
---

# Immunoinformatics Analysis Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for the mathematical formulations, format harmonization rules, HLA HED approximation details, and neoantigen prioritization guardrails.

## Tool Overview

```
AIRR / 10x / TRUST4 / MiXCR   ──imm_repertoire_inspect──▶ harmonized AIRR TSV + QC + spectratype
Harmonized / raw repertoire   ──imm_diversity───────────▶ Shannon/Pielou/D50, homeostasis, V-J usage
Multi-sample repertoire       ──imm_clonotype_tracking──▶ Morisita-Horn/Jaccard, motifs, Fisher clones
OptiType/arcasHLA/seq2HLA/LA  ──imm_hla_summary─────────▶ 2-field consensus, Class I homozygosity, HED
MHC binding / somatic peptide ──imm_neoantigen_prioritize▶ IC50, DAI, anchor/TCR, TPM/VAF gates, rank
```

## Hard Red Lines

1. **Inspect before aggregating**: run `imm_repertoire_inspect` first to verify productive %, out-of-frame/stop-codon rates, and chain composition (`TRA`/`TRB`/`IGH`).
2. **Isolate target loci for diversity & tracking**: compute diversity and longitudinal expansion on a single productive chain (default `TRB`) unless pooling is explicitly justified.
3. **Use two-sided Fisher's exact + BH FDR for clone tracking**: never call longitudinal expansion/contraction from raw fold-change without exact count testing and multiple-testing correction.
4. **Report approximate HED & Class I homozygosity explicitly**: standardize HLA alleles to 2-field (`HLA-X*XX:XX`), flag Class I homozygosity (`HLA-A/B/C`), and label pocket-Grantham HED clearly as an approximation.
5. **Gate & disclose neoantigen predictions**: enforce `ic50_mt_nm <= ic50_cutoff_nm` and `tpm >= min_tpm` gates, separate anchor-driven DAI from TCR-contact foreignness, and state clearly when built-in PSSM estimation is used in place of NetMHCpan/MHCflurry.
