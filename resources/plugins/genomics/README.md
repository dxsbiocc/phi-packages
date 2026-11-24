# genomics

Downstream somatic cancer genomics plugin (`toolPrefix: gen`, private pure-conda
environment `plugin:genomics`).

## Scope & division of work

- **Upstream (out of scope here)**: raw FASTQ alignment, duplicate marking,
  germline/somatic variant calling, and VCF-to-MAF annotation (`nf-core/sarek`,
  Mutect2, Strelka2, VEP / `vcf2maf`, CNVkit / ASCAT).
- **Live database queries**: GDC, cBioPortal, and PDC lookups run via the
  `cancer-omics` MCP connector.
- **This plugin (`gen_*`)**:
  1. `gen_inspect` — validate GDC MAF fields, non-synonymous/silent breakdown,
     VAF/depth coverage, clinical table alignment, and optional CNV `.seg` file.
  2. `gen_summary` — cohort mutation summary, 6-class substitution spectrum and
     Ti/Tv ratio, top mutated genes, MAF summary dashboard, and oncoplot
     (with `Multi_Hit` tracking).
  3. `gen_tmb_ith` — FoCR-compliant TMB (`mutations / capture_size_mb` with
     explicit target territory and VAF/coding filters) and Mroz/Rocco `MATH`
     intratumor heterogeneity score (optionally excluding CNV-altered segments).
  4. `gen_interactions` — pairwise Fisher's exact test for somatic
     co-occurrence/mutual exclusivity (+ BH FDR), Sanchez-Vega (Cell 2018) 10
     canonical oncogenic signaling pathways, and cohort/clinical-group
     differential mutation testing (Fisher + OR + forest plot).
  5. `gen_signatures` — SBS96 trinucleotide spectrum extraction and
     `scipy.optimize.nnls` reference signature refitting (e.g., COSMIC v3.4)
     with cosine similarity evaluation; cleanly falls back to 6-class Ti/Tv
     with an explicit gap note when trinucleotide context is absent.
  6. `gen_cnv` — CBS/CNVkit `.seg` parsing, per-sample Fraction Genome Altered
     (FGA), recurrent gain/loss chromosome/segment summary, and optional
     clinical group comparison.
