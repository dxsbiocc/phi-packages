---
name: guide-design-workflow
description: >-
  CRISPR guide RNA design, off-target specificity scanning, and pooled sgRNA
  library QC: design SpCas9, SpCas9-NG, SaCas9, or Cas12a guides on both strands
  with sequence-context efficiency scoring, Bae et al. (2014) microhomology
  out-of-frame knockout scoring, CBE/ABE base-editor window & bystander/stop
  auditing, and Golden Gate cloning oligos; scan reference FASTAs for 0..4
  mismatch + NAG/NGA off-targets with Hsu-Zhang MIT and Doench 2016 CFD scores;
  and audit sgRNA libraries before oligo synthesis or mageck count. Use for
  guide RNA design, sgRNA, Cas9, Cas12a, base editing, off-target CFD/MIT,
  向导RNA设计, 脱靶评分, 文库质控.
phi:
  environment: plugin:crispr
  attachTo:
    - Crispr
  scripts:
    - name: design_guides
      description: >-
        Design CRISPR guide RNAs from target FASTA sequences across SpCas9,
        SpCas9-NG, SaCas9, or Cas12a nucleases in knockout, CRISPRi, CRISPRa,
        CBE, or ABE mode, computing efficiency score, Bae et al. MMEJ
        out-of-frame KO score, Pol III/restriction enzyme flags, CBE/ABE window
        bystander & stop-codon annotations, and ready-to-order Golden Gate oligos.
      run: [python, ./scripts/crispr_design_guides.py]
      args:
        type: object
        additionalProperties: false
        required: [target_fasta, output_dir]
        properties:
          target_fasta:
            type: string
            format: input-path
            description: Input FASTA containing one or more target transcripts, exons, or genomic loci.
          output_dir:
            type: string
            format: project-path
            description: Directory for designed_guides.tsv and guide_design_overview.png.
          nuclease:
            type: string
            enum: [SpCas9, SpCas9-NG, SaCas9, Cas12a]
            description: Nuclease PAM and spacer architecture (default SpCas9).
          mode:
            type: string
            enum: [knockout, CRISPRi, CRISPRa, CBE, ABE]
            description: Editing modality (default knockout).
          cloning_backbone:
            type: string
            enum: [lentiCRISPRv2_BsmBI, pX330_BbsI, pLKO5_BsmBI, none]
            description: Golden Gate vector overhangs for oligo generation (default lentiCRISPRv2_BsmBI).
          avoid_enzymes:
            type: array
            items:
              type: string
            description: Restriction enzyme names whose recognition sites should be flagged (default BsmBI, BbsI).
          gc_min:
            type: number
            description: Minimum acceptable spacer GC fraction (default 0.35).
          gc_max:
            type: number
            description: Maximum acceptable spacer GC fraction (default 0.75).
          top_n:
            type: integer
            minimum: 1
            description: Number of top-ranked guides highlighted per target in JSON summary (default 10).
      approval: write
      output: ./schemas/design-guides-result.json
      timeoutSeconds: 300
    - name: offtarget_scan
      description: >-
        Perform a PAM-anchored mismatch off-target search (0..4 mismatches,
        optional NAG/NGA non-canonical PAMs) across a reference/panel FASTA,
        computing per-site and aggregate Hsu-Zhang MIT and Doench 2016 CFD
        specificity scores and flagging seed-perfect off-targets.
      run: [python, ./scripts/crispr_offtarget_scan.py]
      args:
        type: object
        additionalProperties: false
        required: [guides, reference_fasta, output_dir]
        properties:
          guides:
            type: string
            format: input-path
            description: TSV/CSV of candidate guides (designed_guides.tsv or 2-column guide_id, spacer).
          reference_fasta:
            type: string
            format: input-path
            description: Reference, locus, transcriptome, or panel FASTA (up to ~100 MB; >250 MB should use nf-core/modules/bowtie or Cas-OFFinder).
          output_dir:
            type: string
            format: project-path
            description: Directory for offtarget_sites.tsv, guide_specificity_summary.tsv, and offtarget_specificity.png.
          nuclease:
            type: string
            enum: [SpCas9, SpCas9-NG, SaCas9, Cas12a]
            description: Nuclease PAM architecture (auto-detected from guides table if omitted, default SpCas9).
          max_mismatches:
            type: integer
            minimum: 0
            maximum: 4
            description: Maximum mismatches allowed in the protospacer (default 3).
          canonical_pam_only:
            type: boolean
            description: Restrict SpCas9 off-target search strictly to canonical NGG PAMs, excluding non-canonical NAG and NGA PAMs (default false).
      approval: write
      output: ./schemas/offtarget-scan-result.json
      timeoutSeconds: 600
    - name: library_check
      description: >-
        Audit a pooled CRISPR sgRNA library (TSV, CSV, or FASTA) prior to oligo
        pool synthesis or mageck count for duplicate IDs/spacers, invalid IUPAC
        bases, length uniformity, guides-per-gene coverage, non-targeting control
        fraction, Pol III TTTT terminators, cloning enzyme collisions, and GC outliers.
      run: [python, ./scripts/crispr_library_check.py]
      args:
        type: object
        additionalProperties: false
        required: [library, output_dir]
        properties:
          library:
            type: string
            format: input-path
            description: sgRNA library TSV/CSV/FASTA table.
          output_dir:
            type: string
            format: project-path
            description: Directory for library_qc_per_gene.tsv, library_flagged_guides.tsv, and library_qc_summary.png.
          sgrna_col:
            type: string
            description: sgRNA ID column name (auto-detected if omitted).
          seq_col:
            type: string
            description: Spacer sequence column name (auto-detected if omitted).
          gene_col:
            type: string
            description: Target gene column name (auto-detected if omitted).
          control_prefix:
            type: string
            description: Comma-separated prefixes/substrings identifying non-targeting or safe-harbor controls (default NonTargeting,SafeHarbor,NTC,CTRL).
          avoid_enzymes:
            type: array
            items:
              type: string
            description: Restriction enzyme recognition motifs to flag inside spacers (default BsmBI, BbsI).
      approval: write
      output: ./schemas/library-check-result.json
      timeoutSeconds: 300
---

# Guide RNA Design, Off-Target Scoring & Library QC Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for nuclease PAM geometries, Bae et al. MMEJ frameshift equations, MIT/CFD off-target scoring, and cloning oligo rules.

## Workflow Chain

```
Target FASTA         ──crispr_design_guides───▶ Candidate guides + MMEJ KO / CBE-ABE window + Golden Gate oligos
Guides + Ref FASTA   ──crispr_offtarget_scan──▶ 0..4 mismatch + NAG/NGA scan → MIT & CFD specificity (0..100)
Pooled sgRNA library ──crispr_library_check───▶ Pre-synthesis / pre-mageck-count audit (dups, TTTT, NTCs, enzymes)
Whole genome (>250MB)──nf-core/modules/bowtie─▶ Genome-wide short-read / oligo alignment wrapper
```

## Hard Red Lines

1. **Knockout target position**: in `knockout` mode, flag guides cutting `<5%` (N-terminal alternative ATG re-initiation risk) or `>85%` (C-terminal NMD escape risk) along the target CDS.
2. **Cloning & Pol III safety**: reject or flag `TTTT` Pol III termination motifs and internal `BsmBI`/`BbsI`/`SapI` cloning enzyme sites before oligo synthesis.
3. **Specificity gate**: treat `MIT_specificity < 50` or `CFD_specificity < 50` or any off-target with `0` PAM-proximal seed mismatches as high off-target risk.
4. **Whole-genome scale**: `crispr_offtarget_scan` is optimized for locus, panel, amplicon, and transcriptome FASTAs (<=100 MB; hard guard at 250 MB). For full mammalian genomes, use `nf-core/modules/bowtie` (`bowtie-build` + `bowtie-align`) or external `Cas-OFFinder`.
