# crispr · 0.2.0

CRISPR gene-editing design (knockout, CRISPRi/a, base editing, Prime Editing), off-target specificity auditing, pooled functional & chemogenetic screening, and post-editing amplicon outcome quantification plugin (`toolPrefix: crispr`, agent `Crispr`) in a locked offline Python environment (`plugin:crispr`: Python 3.13 + `pandas/numpy/scipy/scikit-learn/matplotlib-base`).

## Division of Work: Nextflow Wrappers vs. `crispr_*` Plugin Tools

Phi separates **heavy raw-read / genome-wide processing** (run through vendored Nextflow wrappers in `resources/wrappers/modules/nf-core/`) from **fast, deterministic downstream design, QC, statistical inference, and editing quantification** (run in `plugin:crispr`):

- **Upstream Nextflow Wrappers (`wrapper_search` / `wrapper_plan_run` / `wrapper_run`)**:
  - `nf-core/modules/mageck-count`: map pooled CRISPR screen FASTQ files against an sgRNA library to produce raw sgRNA count tables (`.count.txt`).
  - `nf-core/modules/mageck-test` & `nf-core/modules/mageck-mle`: upstream MAGeCK RRA / MLE CLI alternatives when multi-condition design matrices or official MAGeCK CLI outputs are needed.
  - `nf-core/modules/crisprcleanr-normalize`: CNV-amplification false-positive correction on sgRNA raw counts before hit calling in cancer cell lines with copy-number amplifications.
  - `nf-core/modules/crispresso2`: deep-sequencing amplicon FASTQ demultiplexing, alignment, and allele frequency table generation (`Alleles_frequency_table_around_sgRNA_*.txt`).
  - `nf-core/modules/bowtie` & `nf-core/modules/bowtie2`: whole-genome reference indexing (`bowtie-build` / `bowtie2-build`) and genome-wide short-read/oligo alignment (`bowtie-align` / `bowtie2-align`).

- **Local Deterministic `crispr_*` Tools (`plugin:crispr`, 8 tools across 2 skills)**:
  - **Skill 1 (`guide-design-workflow`, 4 tools)**:
    - `crispr_design_guides`: multi-nuclease (`SpCas9`, `SpCas9-NG`, `SaCas9`, `Cas12a`) guide RNA design on both strands, Doench/sequence-context efficiency scoring, Bae et al. (2014) microhomology out-of-frame KO scoring, `CBE`/`ABE` editing-window + bystander + premature stop auditing, and Golden Gate cloning oligo generation (`lentiCRISPRv2_BsmBI`, `pX330_BbsI`, `pLKO5_BsmBI`).
    - `crispr_design_pegrna`: Prime Editing (`PE2` / `PE3` / `PE3b`) pegRNA, epegRNA (`3' tevopreQ1` motif), and secondary nicking sgRNA (`ngRNA`) design from wild-type vs. edited target sequences, featuring SantaLucia PBS $T_m$ scoring, RTT `5'-C` Cas9-scaffold pairing guard, synonymous PAM-disrupting codon discovery, and PE3b edit-specific nickase prioritization.
    - `crispr_offtarget_scan`: PAM-anchored mismatch off-target search (`0..4` mismatches + non-canonical `NAG`/`NGA` PAMs) against reference/panel FASTAs with Hsu-Zhang **MIT** and Doench 2016 **CFD** specificity scores.
    - `crispr_library_check`: pre-synthesis / pre-`mageck-count` audit of pooled sgRNA libraries (duplicate IDs/spacers, invalid IUPAC bases, length uniformity, guides-per-gene coverage, NTC control representation, `TTTT` terminators, restriction site collisions).
  - **Skill 2 (`screen-editing-workflow`, 4 tools)**:
    - `crispr_screen_qc`: pooled screen count QC (median reads/sgRNA, zero-dropout %, Lorenz curve & Gini index, DESeq median-of-ratios size factors, sample Pearson correlation & PCA, and essential vs non-essential **NNMD** + ROC-AUC).
    - `crispr_screen_hits`: normalized sgRNA-level $\log_2\text{FC}$ + moderated variance modeling, gene-level $\alpha$-RRA positive/negative selection scoring with seeded empirical permutation calibration and Benjamini–Hochberg FDR.
    - `crispr_drugz`: Collimore & Hart (2019) **DrugZ** chemogenetic drug-vs-vehicle pooled CRISPR screen analysis (10M RPM normalization, abundance-ordered empirical-Bayes monotonic sliding-window variance estimation, paired/unpaired sgRNA $z$-scores, gene-level `sumZ`/`normZ`, and two-sided BH FDR for synthetic-lethal sensitizers & drug-resistance suppressors).
    - `crispr_editing_quant`: post-editing outcome quantification from raw amplicon allele tables/FASTAs or `crispresso2` allele frequency tables — separating **Unmodified**, **Frameshift Indel** (productive KO), **In-frame Indel**, **Pure Substitution**, **HDR**, and **Imperfect HDR**, plus indel spectrum and protospacer base-conversion (`CBE`/`ABE` target vs bystander) matrices.
