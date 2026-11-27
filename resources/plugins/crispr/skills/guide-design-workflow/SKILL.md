---
name: guide-design-workflow
description: >-
  CRISPR guide RNA & Prime Editing pegRNA design, off-target specificity
  scanning, and pooled sgRNA library QC: design SpCas9, SpCas9-NG, SaCas9, or
  Cas12a guides on both strands with sequence-context efficiency scoring, Bae et
  al. (2014) microhomology out-of-frame knockout scoring, CBE/ABE base-editor
  window & bystander/stop auditing, and Golden Gate cloning oligos; design
  PE2/PE3/PE3b pegRNAs, epegRNAs (tevopreQ1), and nicking ngRNAs; scan reference
  FASTAs for 0..4 mismatch + NAG/NGA off-targets with Hsu-Zhang MIT and Doench
  2016 CFD scores; and audit sgRNA libraries before oligo synthesis or mageck
  count. Use for guide RNA design, sgRNA, pegRNA, epegRNA, Prime Editing, Cas9,
  Cas12a, base editing, off-target CFD/MIT, 向导RNA设计, 引导编辑, 脱靶评分, 文库质控.
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
    - name: design_pegrna
      description: >-
        Design Prime Editing (PE2, PE3, and PE3b) pegRNAs, optional 3'-tevopreQ1
        epegRNAs, and opposite-strand nicking sgRNAs (ngRNAs) from wild-type and
        edited DNA sequences (>=60 bp), scoring SantaLucia PBS Tm, RTT 5'-C
        scaffold-pairing penalties, Pol III TTTT motifs, and synonymous
        PAM-disrupting codon hints.
      run: [python, ./scripts/crispr_design_pegrna.py]
      args:
        type: object
        additionalProperties: false
        required: [wt_seq, edited_seq, output_dir]
        properties:
          wt_seq:
            type: string
            description: Unambiguous wild-type DNA sequence (>= 60 bp) flanking the desired edit (5'->3').
          edited_seq:
            type: string
            description: Desired post-edit DNA sequence (>= 60 bp) sharing 5' and 3' flanks with wt_seq (5'->3').
          output_dir:
            type: string
            format: project-path
            description: Directory for pegrna_candidates.tsv, ngrna_candidates.tsv, and pegrna_architecture.png.
          nuclease:
            type: string
            enum: [SpCas9, SpCas9-NG]
            description: Prime editor nickase PAM architecture (default SpCas9).
          pbs_min:
            type: integer
            minimum: 7
            maximum: 20
            description: Minimum Primer Binding Site (PBS) length in nt (default 10).
          pbs_max:
            type: integer
            minimum: 8
            maximum: 22
            description: Maximum Primer Binding Site (PBS) length in nt (default 15).
          rtt_homology_min:
            type: integer
            minimum: 5
            maximum: 30
            description: Minimum post-edit 3' flap homology length in the RT template in nt (default 8).
          rtt_max:
            type: integer
            minimum: 10
            maximum: 50
            description: Maximum total RT template (RTT) length in nt (default 32).
          add_epgrna_motif:
            type: boolean
            description: Append the structured 3' tevopreQ1 pseudoknot motif (with 4-nt linker) to protect epegRNAs from 3' exonucleases (default false).
          cds_frame_offset:
            type: integer
            enum: [0, 1, 2]
            description: Optional 0-based reading-frame start offset on edited_seq to suggest in-frame synonymous PAM-disrupting mutations.
      approval: write
      output: ./schemas/design-pegrna-result.json
      timeoutSeconds: 300
---

# Guide RNA & pegRNA Design, Off-Target Scoring & Library QC Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for nuclease PAM geometries, Bae et al. MMEJ frameshift equations, Prime Editing PBS/RTT/PE3b rules, MIT/CFD off-target scoring, and cloning oligo rules.

## Workflow Chain

```
Target FASTA         ──crispr_design_guides───▶ Candidate guides + MMEJ KO / CBE-ABE window + Golden Gate oligos
WT + Edited DNA      ──crispr_design_pegrna───▶ PE2/PE3/PE3b pegRNAs + epegRNA (tevopreQ1) + ngRNAs + PBS Tm
Guides + Ref FASTA   ──crispr_offtarget_scan──▶ 0..4 mismatch + NAG/NGA scan → MIT & CFD specificity (0..100)
Pooled sgRNA library ──crispr_library_check───▶ Pre-synthesis / pre-mageck-count audit (dups, TTTT, NTCs, enzymes)
Whole genome (>250MB)──nf-core/modules/bowtie─▶ Genome-wide short-read / oligo alignment wrapper
```

## Hard Red Lines

1. **Knockout target position**: in `knockout` mode, flag guides cutting `<5%` (N-terminal alternative ATG re-initiation risk) or `>85%` (C-terminal NMD escape risk) along the target CDS.
2. **Prime Editing RTT & PE3b rules**: avoid pegRNA RTTs starting with `5'-C` (pairs with SpCas9 scaffold G81), target PBS $T_m \approx 37\text{--}42^\circ\text{C}$, and prefer edit-specific `PE3b` nicking guides over distal `PE3` nicking guides to minimize DSB indels.
3. **Cloning & Pol III safety**: reject or flag `TTTT` Pol III termination motifs and internal `BsmBI`/`BbsI`/`SapI` cloning enzyme sites before oligo synthesis.
4. **Specificity gate**: treat `MIT_specificity < 50` or `CFD_specificity < 50` or any off-target with `0` PAM-proximal seed mismatches as high off-target risk.
5. **Whole-genome scale**: `crispr_offtarget_scan` is optimized for locus, panel, amplicon, and transcriptome FASTAs (<=100 MB; hard guard at 250 MB). For full mammalian genomes, use `nf-core/modules/bowtie` (`bowtie-build` + `bowtie-align`) or external `Cas-OFFinder`.
