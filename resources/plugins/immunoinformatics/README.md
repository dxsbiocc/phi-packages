# immunoinformatics · 0.1.0

Downstream TCR/BCR immune repertoire, HLA genotyping consensus, and neoantigen prioritization plugin (`toolPrefix: imm`, agent `Immunoinformatics`) running in a locked MIT/BSD `conda-forge` environment (`plugin:imm`: Python 3.13 + `pandas`/`numpy`/`scipy`/`scikit-learn`/`matplotlib-base`).

All code in `plugin:imm` is 100% original MIT-licensed Python with zero AGPL (`immunarch`/`immcantation`) or NPOSL (`pVACtools`) dependencies.

## Division of Work: Upstream `nf-core` Wrappers vs. `imm_*` Tools

| Stage / Task | Heavy Upstream Wrapper (`nf-core/modules/*`) | Deterministic Downstream Tool (`plugin:imm`) |
| --- | --- | --- |
| Raw RNA-seq / scRNA-seq TCR & BCR de novo assembly | `nf-core/modules/trust4` | — |
| Class I HLA typing from DNA/RNA BAMs (ILP) | `nf-core/modules/optitype` | — |
| Class I & II HLA typing from RNA-seq BAMs | `nf-core/modules/arcashla` | — |
| Class I & II HLA typing & expression from RNA-seq FASTQs | `nf-core/modules/seq2hla` | — |
| Graph-based Class I & II HLA typing from WGS/WES/long reads | `nf-core/modules/hlala` | — |
| AIRR / 10x Cell Ranger V(D)J / TRUST4 / MiXCR harmonization & QC | — | `imm_repertoire_inspect` |
| Shannon/Pielou/Simpson/D50 diversity, homeostasis tiers & V-J usage | — | `imm_diversity` |
| Morisita-Horn/Jaccard overlap, 4-mer motifs & Fisher clonotype tracking | — | `imm_clonotype_tracking` |
| Multi-tool 2-field HLA consensus, LOH/homozygosity & approximate HED | — | `imm_hla_summary` |
| MHC IC50, agretopicity (DAI), anchor/foreignness & TESLA prioritization | — | `imm_neoantigen_prioritize` |

## Tools (`imm_*`)

- `imm_repertoire_inspect`: Harmonizes AIRR, 10x `filtered_contig_annotations.csv`, TRUST4 reports, and MiXCR clone tables into a standardized AIRR-compatible TSV (`harmonized_repertoire.tsv`), per-sample QC summary (productive %, out-of-frame/stop %, chain breakdown, CDR3 length, Kyte-Doolittle hydropathy), and spectratype figure.
- `imm_diversity`: Computes Shannon entropy (bits/nats), Pielou clonality (`1 - H / ln(S)`), Gini-Simpson, Inverse Simpson, `D50`, clonal homeostasis space (`Hyperexpanded >=1%`, `Large 0.1-1%`, `Medium 0.01-0.1%`, `Small <0.01%`), V-J gene pairing usage matrix, and optional Mann-Whitney / Kruskal-Wallis + BH FDR group tests.
- `imm_clonotype_tracking`: Computes pairwise Morisita-Horn and Jaccard CDR3AA+V overlap matrices, public/shared clonotypes and enriched 4-mer CDR3 motifs, plus two-sided Fisher's exact test + BH FDR calling significantly `expanded` vs `contracted` clones between baseline and follow-up samples.
- `imm_hla_summary`: Standardizes OptiType, arcasHLA, seq2HLA, HLA-LA, and cohort HLA tables to 2-field (`4-digit`, e.g. `HLA-A*02:01`) across Class I (`HLA-A/B/C`) and Class II (`HLA-DRB1/DQB1/DPB1`), flags Class I homozygosity, and computes approximate pocket-Grantham HLA Evolutionary Divergence (HED).
- `imm_neoantigen_prioritize`: Ingests upstream MHC binding tables or somatic peptide tables (with built-in anchor-motif PSSM IC50 estimation when `ic50_mt` is absent), computes Differential Agretopicity Index (`DAI = IC50_wt / IC50_mt`), anchor vs TCR-contact mutation status, foreignness/hydropathy, expression (`TPM`) & clonality (`VAF`/`CCF`) gates, and composite `immunogenicity_score` (`0..100`).
