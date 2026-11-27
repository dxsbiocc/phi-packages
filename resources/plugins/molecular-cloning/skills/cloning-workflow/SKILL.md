---
name: cloning-workflow
description: >-
  End-to-end molecular cloning, primer design, codon optimization, construct
  assembly, and Sanger verification workflow: design PCR/qPCR/Gibson/GoldenGate/SDM
  primers with SantaLucia (1998) + Owczarzy (2008) Tm and hairpin/dimer dG37
  scoring; optimize coding sequences for E. coli, HEK293, CHO, or yeast while
  eliminating restriction sites; assemble constructs via Golden Gate, Gibson, or
  Restriction Ligation with virtual gel verification; and verify clones from
  binary ABIF .ab1, FASTQ, or FASTA Sanger traces including mixed-peak
  deconvolution. Use for primer design, Gibson assembly, Golden Gate cloning,
  codon optimization, plasmid map, Sanger .ab1 verification, 引物设计, 密码子优化,
  载体构建, 一代测序验证.
phi:
  environment: plugin:mc
  attachTo:
    - MolecularCloning
  scripts:
    - name: design_primers
      description: >-
        Design primer pairs for PCR, qPCR, Gibson, GoldenGate, or
        SiteDirectedMutagenesis from FASTA or GenBank templates using
        SantaLucia (1998) nearest-neighbor thermodynamics and Owczarzy et al.
        (2008) Mg2+/dNTP/Na+ salt correction, scoring GC%, 3' clamp,
        homopolymers, hairpin/self-dimer/cross-dimer dG37, and optional
        background 3'-12mer mispriming.
      run: [python, ./scripts/mc_design_primers.py]
      args:
        type: object
        additionalProperties: false
        required: [template, output_dir]
        properties:
          template:
            type: string
            format: input-path
            description: Input template DNA in FASTA or GenBank (.gb/.gbk) format.
          output_dir:
            type: string
            format: project-path
            description: Output directory for primer_candidates.tsv and primer_design_map.png.
          mode:
            type: string
            enum: [PCR, qPCR, Gibson, GoldenGate, SiteDirectedMutagenesis]
            description: Primer design modality (default PCR).
          target_start:
            type: integer
            minimum: 1
            description: 1-based start coordinate of target region on template.
          target_end:
            type: integer
            minimum: 1
            description: 1-based end coordinate of target region on template.
          amplicon_min:
            type: integer
            minimum: 30
            description: Minimum amplicon length in bp (default 80 for qPCR, 200 otherwise).
          amplicon_max:
            type: integer
            minimum: 40
            description: Maximum amplicon length in bp (default 180 for qPCR, 1200 otherwise).
          tm_opt:
            type: number
            description: Target annealing Tm in Celsius (default 60.0).
          na_mm:
            type: number
            description: Monovalent Na+/K+ concentration in mM (default 50.0).
          mg_mm:
            type: number
            description: Divalent Mg2+ concentration in mM (default 1.5).
          dntp_mm:
            type: number
            description: Total dNTP concentration in mM (default 0.2).
          primer_nm:
            type: number
            description: Total primer concentration in nM (default 250.0).
          left_homology:
            type: string
            description: 5' Gibson vector overlap homology tail added to forward primer.
          right_homology:
            type: string
            description: 3' Gibson vector overlap homology tail (revcomp added to reverse primer).
          goldengate_enzyme:
            type: string
            enum: [BsaI, BsmBI, BbsI, SapI]
            description: Type IIS enzyme for GoldenGate primer adapters (default BsaI).
          overhang_5p:
            type: string
            description: 4-nt (or 3-nt for SapI) 5' Golden Gate overhang sequence.
          overhang_3p:
            type: string
            description: 4-nt (or 3-nt for SapI) 3' Golden Gate overhang sequence.
          mutation_pos:
            type: integer
            minimum: 1
            description: 1-based mutation position on template for SiteDirectedMutagenesis.
          mutation_seq:
            type: string
            description: Replacement DNA bases at mutation_pos for SiteDirectedMutagenesis.
          background_fasta:
            type: string
            format: input-path
            description: Optional off-target FASTA file to scan for 3'-12mer mispriming hits.
      approval: write
      output: ./schemas/design-primers-result.json
      timeoutSeconds: 300
    - name: codon_optimize
      description: >-
        Optimize a protein or DNA CDS FASTA for host expression (E_coli_K12,
        HEK293_Human, CHO_Cricetulus, S_cerevisiae), maximizing CAI while
        preserving 100% amino-acid identity, smoothing 30-bp GC windows,
        relaxing 5' secondary structure, and eliminating forbidden restriction
        sites, homopolymers, and cryptic E. coli Shine-Dalgarno motifs.
      run: [python, ./scripts/mc_codon_optimize.py]
      args:
        type: object
        additionalProperties: false
        required: [sequence_input, output_dir]
        properties:
          sequence_input:
            type: string
            format: input-path
            description: Input protein FASTA or DNA CDS FASTA.
          output_dir:
            type: string
            format: project-path
            description: Output directory for optimized_cds.fasta, optimized_construct.gb, codon_audit.tsv, and codon_optimization_profile.png.
          is_dna:
            type: boolean
            description: Treat input strictly as a DNA CDS rather than auto-detecting protein vs DNA (default false).
          host:
            type: string
            enum: [E_coli_K12, HEK293_Human, CHO_Cricetulus, S_cerevisiae]
            description: Target expression host codon table (default E_coli_K12).
          avoid_enzymes:
            type: array
            items:
              type: string
            description: Restriction enzyme recognition sites to eliminate (default BsaI, BsmBI, BbsI, EcoRI, BamHI, HindIII, NdeI, XhoI).
          gc_window_min:
            type: number
            description: Minimum 30-bp sliding-window GC fraction (default 0.30).
          gc_window_max:
            type: number
            description: Maximum 30-bp sliding-window GC fraction (default 0.72).
          add_kozak:
            type: boolean
            description: Prepend GCCACC Kozak sequence 5' of start codon in output construct (default false).
          seed:
            type: integer
            description: Deterministic random seed for optimization tie-breaking (default 42).
      approval: write
      output: ./schemas/codon-optimize-result.json
      timeoutSeconds: 300
    - name: assemble_construct
      description: >-
        Virtually assemble DNA fragments and optional vector backbone using
        GoldenGate (Type IIS overhang graph + fidelity & retained-site audit),
        Gibson (15-60 bp terminal overlap stitching & Tm), or
        RestrictionLigation (35+ Type II enzymes with compatible cohesive ends),
        exporting FASTA, annotated GenBank .gb, junction audit TSV, and a
        2-panel plasmid map + virtual 1% agarose gel.
      run: [python, ./scripts/mc_assemble_construct.py]
      args:
        type: object
        additionalProperties: false
        required: [fragments, output_dir, method]
        properties:
          fragments:
            type: string
            format: input-path
            description: Multi-record FASTA or GenBank (.gb) containing ordered DNA fragments.
          output_dir:
            type: string
            format: project-path
            description: Output directory for assembled_construct.fasta, assembled_construct.gb, assembly_junctions.tsv, and construct_map_and_gel.png.
          method:
            type: string
            enum: [GoldenGate, Gibson, RestrictionLigation]
            description: Cloning assembly strategy.
          vector:
            type: string
            format: input-path
            description: Optional vector backbone FASTA or GenBank (.gb) file.
          enzyme:
            type: string
            description: Enzyme name(s) (e.g. BsaI for GoldenGate, or EcoRI,BamHI for RestrictionLigation).
          linear_product:
            type: boolean
            description: Assemble a linear DNA product instead of a circular plasmid (default false).
          construct_name:
            type: string
            description: Name for assembled plasmid/construct in FASTA and GenBank LOCUS (default pAssembled).
      approval: write
      output: ./schemas/assemble-construct-result.json
      timeoutSeconds: 300
    - name: sanger_verify
      description: >-
        Verify assembled constructs or amplicons from binary ABIF .ab1
        chromatograms, FASTQ, or FASTA Sanger reads using Mott Q-score trimming,
        strand auto-orientation with circular origin wrap-around alignment,
        SNV/indel/frameshift calling, and secondary fluorescence peak ratio
        deconvolution.
      run: [python, ./scripts/mc_sanger_verify.py]
      args:
        type: object
        additionalProperties: false
        required: [reference, traces, output_dir]
        properties:
          reference:
            type: string
            format: input-path
            description: Expected construct/amplicon sequence in FASTA or GenBank (.gb/.gbk) format.
          traces:
            type: array
            items:
              type: string
              format: input-path
            description: Paths to Sanger sequencing traces (.ab1, .fastq, or .fasta).
          output_dir:
            type: string
            format: project-path
            description: Output directory for sanger_verification_summary.tsv, sanger_variants.tsv, and sanger_alignment_chromatogram.png.
          trim_error_prob:
            type: number
            description: Mott algorithm error probability cutoff for quality trimming (default 0.05).
          min_secondary_peak_ratio:
            type: number
            description: Minimum I2/(I1+I2) secondary channel ratio in .ab1 traces to call mixed/heterozygous peaks (default 0.15).
      approval: write
      output: ./schemas/sanger-verify-result.json
      timeoutSeconds: 300
---

# Molecular Cloning, Primer Design & Sanger Verification Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for SantaLucia (1998) + Owczarzy (2008) $T_m$ equations, Golden Gate overhang fidelity rules, codon optimization constraints, and ABIF `.ab1` secondary peak thresholds.

## End-to-End Cloning Pipeline

```
Protein / CDS FASTA  ──mc_codon_optimize─────▶ Host-optimized CDS + .gb (0 restriction sites)
Template FASTA / .gb ──mc_design_primers─────▶ PCR/qPCR/Gibson/GoldenGate/SDM primers + Tm/ΔG37
Ordered Parts / Vec  ──mc_assemble_construct─▶ Assembled .gb + junction fidelity + virtual gel
Sanger .ab1 / .fastq ──mc_sanger_verify──────▶ Circular alignment + SNV/indel + mixed-peak check
```

## Hard Red Lines

1. **Annealing vs. full-oligo $T_m$**: For `Gibson` and `GoldenGate` primers, evaluate $|\Delta T_m| \le 2.5^\circ\text{C}$ on the **3' template-annealing core** separately from `full_tm`.
2. **Strict translation invariant**: `mc_codon_optimize` must preserve 100% identical amino-acid translation while clearing all forbidden restriction enzymes and Pol III terminators (`TTTT`).
3. **Golden Gate overhang hygiene**: Flag palindromic overhangs (`ATAT`, `GCGC`, `ACGT`, `TGCA`), duplicate/reverse-complement overhang collisions, 1-bp mismatch crosstalk, and any retained internal Type IIS sites.
4. **Sanger mixed-colony vigilance**: In `.ab1` traces, inspect secondary peak ratios $I_2/(I_1+I_2) \ge 0.15$ alongside SNVs/indels so mixed colonies or heterozygous edits are never misclassified as clean clones.
