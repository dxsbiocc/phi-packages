# Changelog — genomics

## 0.3.0 (2026-11-26)

- Add `gen_hrd_msi`: computes scarHRD/Marquard genomic scar Homologous Recombination Deficiency (`HRD-LOH`, `TAI`, `LST`, and `hrd_sum`) on hg38/hg19 autosomes from allele-specific CNV segments (ASCAT/Sequenza/FACETS/PURPLE), integrates MSIsensor-pro/Mantis `%` unstable sites or somatic MAF frameshift-indel/indel-to-SNV surrogates, and classifies therapeutic biomarker strata (`HRD_pos_MSS`, `MSI_H`, `HR_proficient_MSS`).
- Add `gen_clonality`: estimates per-variant Cancer Cell Fraction (`CCF`) with Binomial multiplicity MLE and Clopper-Pearson 95% confidence intervals under sample purity and local copy number, separates clonal vs subclonal driver mutations, and evaluates Williams et al. (2016) neutral tail linearity ($M(f) \propto 1/f$, $0.12 \le \text{VAF} \le 0.24$).
- Update `Genomics.md`, `genomics-workflow` `SKILL.md`, `analysis-playbook.md`, JSON schemas, and `tests/test_genomics_plugin.py` for 0.3.0 (9 tools).

## 0.2.0 (2026-11-25)

- Add `gen_vcf2maf`: deterministic converter from VEP (`CSQ`), SnpEff (`ANN`), or ANNOVAR annotated somatic VCFs (`.vcf`/`.vcf.gz`, e.g. `nf-core/sarek` outputs) to GDC-compliant MAF tables, including canonical/severity transcript prioritization, Mutect2/Strelka2 VAF & depth extraction, and optional reference FASTA 3-mer `Context` extraction for downstream `gen_signatures` SBS96 refitting.
- Add contract-1.0.0 `.phi-artifact.json` descriptors and `artifacts` lists across all file-producing tools (`gen_vcf2maf`, `gen_summary`, `gen_tmb_ith`, `gen_interactions`, `gen_signatures`, `gen_cnv`) and update output JSON schemas.

## 0.1.0 (2026-11-24)

Initial release: pure conda-forge downstream somatic cancer genomics plugin
(`toolPrefix: gen`).

- `gen_inspect`: GDC MAF validation, non-synonymous/silent breakdown, VAF/depth
  QC, clinical alignment, and optional CNV `.seg` inspection.
- `gen_summary`: cohort variant summary, 6-class substitution spectrum & Ti/Tv,
  MAF summary dashboard, and waterfall oncoplot (`Multi_Hit`).
- `gen_tmb_ith`: FoCR-compliant TMB (explicit `capture_size_mb` + `min_vaf`)
  and Mroz/Rocco `MATH` intratumor heterogeneity score with optional CNV segment
  exclusion and clinical group testing.
- `gen_interactions`: pairwise Fisher co-occurrence/mutual exclusivity (+ BH
  FDR), Sanchez-Vega (Cell 2018) 10 canonical oncogenic pathways, and
  cohort/clinical-group differential mutation testing (Fisher + OR + forest
  plot).
- `gen_signatures`: SBS96 trinucleotide context spectrum and `scipy.optimize.nnls`
  reference signature refitting with cosine similarity diagnostics.
- `gen_cnv`: CBS/CNVkit `.seg` Fraction Genome Altered (FGA), recurrent gain/loss
  chromosome summaries, and clinical group comparison.
