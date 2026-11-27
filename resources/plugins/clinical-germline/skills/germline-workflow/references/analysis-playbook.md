# Clinical Germline & PGx Analysis Playbook

## 0. Division of Work with Upstream `nf-core` Wrappers & Connectors

- **Upstream calling & annotation** runs via `nf-core` wrappers:
  - `nf-core/sarek`: end-to-end germline WGS/WES alignment, BQSR, and joint variant calling.
  - `nf-core/modules/deepvariant` & `nf-core/modules/clair3`: short-read and long-read/hybrid SNV/indel calling.
  - `nf-core/modules/manta`: paired-end / split-read germline structural variant (DEL/DUP/INV/BND) calling.
  - `nf-core/modules/ensemblvep` & `nf-core/modules/snpeff`: transcript consequence (`CSQ`/`ANN`), population AF (`gnomAD_AF`), ClinVar, and in-silico deleteriousness (`REVEL`, `CADD`, `SpliceAI`) annotation.
- **Ontology lookup** uses the `ols4` MCP connector (EMBL-EBI OLS4) to resolve unfamiliar HPO identifiers and parent/child disease phenotypes.
- **Downstream clinical interpretation** runs deterministically in `plugin:germ` (`germ_*`), implemented in 100% original MIT-licensed Python with zero non-commercial `InterVar` or `ANNOVAR` source code.

## 1. Mendelian Inheritance & Phasing Filter (`germ_vcf_filter`)

- **QC & Coding/Splice Gate**:
  - Filters proband carrier genotypes (`0/1`, `1/1`, or hemizygous `1`) by read depth (`DP >= min_dp`, default `10`) and genotype quality (`GQ >= min_gq`, default `20`) when those FORMAT/tabular fields are present.
  - Retains protein-altering or canonical/region splice consequences (`stop_gained`, `frameshift`, `splice_acceptor`, `splice_donor`, `splice_region`, `missense`, `inframe_insertion`, `inframe_deletion`, `start_lost`, `stop_lost`, `protein_altering`), variants with `spliceai_max >= 0.20`, or ClinVar Pathogenic/Likely Pathogenic entries.
- **Mendelian Fit**:
  - `de_novo`: Proband heterozygous (`0/1`), both parents homozygous reference (`0/0`) with parental QC verified, and `gnomad_af <= max_af_dominant` (default `0.0001`).
  - `autosomal_recessive_hom`: Proband homozygous alternate (`1/1`) with `gnomad_af <= max_af_recessive` (default `0.01`), confirmed biparental carrier (`0/1` + `0/1`) in trios when parental genotypes are present.
  - `compound_heterozygous`: `>= 2` rare (`<= max_af_recessive`) coding/splice variants in the **same gene** in the proband. In trio mode, requires one variant exclusively maternal and the other exclusively paternal (`phased_in_trans`; variants sharing the same parent are flagged `in_cis_excluded`). In singleton mode without parental genotypes, flagged **`unphased_candidate_comp_het`**.
  - `x_linked`: `chrX` variants evaluated against `proband_sex` (hemizygous/homozygous in `male`, heterozygous/homozygous in `female`).
  - `autosomal_dominant`: Rare heterozygous variant with `gnomad_af <= max_af_dominant`.

## 2. ACMG/AMP 2015 + Tavtigian Bayesian Classification (`germ_acmg_prioritize`)

Implements **Richards et al. (2015)** criteria calibrated with **Pejaver et al. (2022)** ClinGen computational thresholds and **Tavtigian et al. (2018, 2020)** Bayesian point scaling ($P_0 = 0.10$, $O_{\text{PVS}} = 350$):

$$\text{Odds}_{\text{post}} = \frac{P_0}{1 - P_0} \cdot O_{\text{PVS}}^{\text{points}/8}, \qquad P_{\text{post}} = \frac{\text{Odds}_{\text{post}}}{1 + \text{Odds}_{\text{post}}}$$

- **Pathogenic Criteria**:
  - `PVS1` (`+8`): Null variant (`stop_gained`, `frameshift`, `splice_acceptor`, `splice_donor`, `start_lost`) in a LoF-intolerant gene (`pLI >= 0.9`, `loeuf <= 0.35`, or curated/user-supplied haploinsufficient gene list).
  - `PS5_ClinVar` (`+4`): Established ClinVar `Pathogenic` assertion without conflicting interpretations or benign assertions.
  - `PS2` (`+4`): Trio-confirmed `de_novo` inheritance.
  - `PM2_Supporting` (`+1`): Absent or ultra-rare in gnomAD (`gnomad_af == 0` or `< 1e-5` dominant / `< 1e-4` recessive).
  - `PM3` (`+2`): Phased *in trans* `compound_heterozygous` in a recessive gene.
  - `PM4` (`+2`): In-frame insertion/deletion in a non-repeat coding region.
  - `PP3` (Pejaver 2022, non-redundant max tier): `PP3_Strong` (`+4`: `REVEL >= 0.932` or `SpliceAI >= 0.80` or `CADD >= 29.8`), `PP3_Moderate` (`+2`: `REVEL >= 0.773` or `SpliceAI >= 0.50` or `CADD >= 28.1`), `PP3_Supporting` (`+1`: `REVEL >= 0.644` or `SpliceAI >= 0.20` or `CADD >= 25.3`).
- **Benign Criteria**:
  - `BA1` (`-8`, stand-alone benign): `gnomad_af > 0.05`.
  - `BS1` (`-4`): `gnomad_af > 0.01`.
  - `BP4` (when `SpliceAI < 0.10` or absent): `BP4_Moderate` (`-2`: `REVEL <= 0.183`), `BP4_Supporting` (`-1`: `REVEL <= 0.290` or `CADD <= 15.0`).
  - `BP7` (`-1`): Synonymous variant with `SpliceAI < 0.10` (or absent).
- **5-Tier Classification**: `Pathogenic` (`>= 10`), `Likely_pathogenic` (`6..9`), `VUS` (`0..5`), `Likely_benign` (`-6..-1`), `Benign` (`<= -7` or `BA1`).
- **ACMG SF v3.2**: Independently flags Pathogenic/Likely Pathogenic variants in the 81 ACMG SF v3.2 medically actionable genes (`acmg_sf_reportable`).

## 3. IC-Weighted HPO Phenotype Ranking (`germ_phenotype_rank`)

- Resolves user-supplied `HP:XXXXXXX` IDs and English/Chinese clinical keywords against a built-in ~60-gene multi-system rare-disease HPO knowledge base (epilepsy/neurodevelopmental, cardiomyopathy/arrhythmia, connective/skeletal, metabolic/lysosomal, immunodeficiency, renal/ciliopathy) or custom `genes_to_phenotype.txt`.
- Weights each HPO term $t$ by corpus Information Content:
  $$\text{IC}(t) = -\ln\left(\frac{n_t + 0.5}{N_{\text{genes}} + 1}\right)$$
  computing `phenotype_sim` $\in [0, 1]$ from IC-weighted recall and Jaccard overlap (with semantic domain credit for related terms).
- Combines variant pathogenicity (`0..50`), `phenotype_sim` (`0..40`), and Mendelian inheritance fit (`0..10`) into an Exomiser-inspired `combined_diagnostic_score` (`0..100`).

## 4. Clinical Pharmacogenomics (`germ_pgx_star_allele`)

- Resolves star-allele diplotypes directly or from defining rsIDs/genotypes across 9 CPIC Level A pharmacogenes: `CYP2C19` (`*2`, `*3`, `*17`), `CYP2D6` (`*4`, `*6`, `*10`, `*41`), `CYP2C9` (`*2`, `*3`), `VKORC1` (`rs9923231`), `SLCO1B1` (`*5`), `DPYD` (`*2A`, `*13`, `c.2846A>T`, `HapB3`), `TPMT` (`*2`, `*3A`, `*3B`, `*3C`), `NUDT15` (`*3`), and `UGT1A1` (`*80/*28`, `*6`).
- Computes CPIC **Activity Scores**, standardized phenotypes (`Ultrarapid_Metabolizer`, `Rapid_Metabolizer`, `Normal_Metabolizer`, `Intermediate_Metabolizer`, `Poor_Metabolizer`, `Decreased_Function`), and structured drug alerts (`contraindicated_or_avoid`, `dose_reduction_required`, `alternative_therapy`, `standard_dosing`).

## 5. Germline CNV/SV + SNV Integration (`germ_cnv_sv_integrate`)

- Evaluates `DEL`, `DUP`, `INV`, and `BND` calls under **ClinGen/ACMG 2020 (Riggs et al.)** dosage-sensitivity criteria (`HI` haploinsufficiency and `TS` triplosensitivity scores, exon/coding disruption, and copy number).
- When `--snv_candidates` is supplied, cross-references rare coding/splice SNVs against heterozygous gene/exon deletions in the same gene to uncover *trans*-acting **SNV + exon-DEL** compound heterozygotes in recessive disorders (`snv_cnv_compound_hits.tsv`).
