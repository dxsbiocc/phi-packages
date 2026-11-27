# genomics

Downstream somatic cancer genomics plugin (`0.3.0`, `toolPrefix: gen`, private
pure-conda environment `plugin:genomics`, 9 tools).

## Scope & division of work

- **Upstream (out of scope here)**: raw FASTQ alignment, duplicate marking,
  germline/somatic variant calling, VCF-to-MAF annotation, and MSI site scanning
  (`nf-core/sarek`, Mutect2, Strelka2, VEP / `vcf2maf`, CNVkit / ASCAT, MSIsensor-pro).
- **Live database queries**: GDC, cBioPortal, and PDC lookups run via the
  `cancer-omics` MCP connector.
- **This plugin (`gen_*`, 9 tools)**:
  0. `gen_vcf2maf` — convert annotated somatic VCFs (`.vcf`/`.vcf.gz` from
     `nf-core/sarek` VEP `CSQ`, SnpEff `ANN`, or ANNOVAR) to GDC-compliant MAF
     tables with severity-prioritized transcript selection, VAF/depth extraction,
     and optional FASTA 3-mer `Context` extraction for SBS96 signatures.
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
  7. `gen_hrd_msi` — scarHRD/Marquard genomic scar HRD scoring (`HRD-LOH`,
     `TAI`, `LST`, `hrd_sum` against clinical cutoff `42`) from allele-specific
     CNV segments (ASCAT/Sequenza/FACETS/PURPLE) combined with MSIsensor-pro /
     Mantis MSI tables or MAF indel/frameshift surrogates for biomarker
     stratification (`HRD_pos_MSS`, `MSI_H`, `HR_proficient_MSS`).
  8. `gen_clonality` — McGranahan/PyClone-inspired Cancer Cell Fraction (`CCF`)
     estimation with Binomial multiplicity MLE, Clopper-Pearson 95% CI,
     clonal vs subclonal driver attribution, and Williams et al. (2016)
     neutral subclonal tail ($M(f) \propto 1/f$) $R^2$ testing.
