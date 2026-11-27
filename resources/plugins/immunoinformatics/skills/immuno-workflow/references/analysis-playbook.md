# Immunoinformatics Analysis Playbook

## 0. Division of Work with Upstream `nf-core` Wrappers

- **Upstream alignment / assembly / raw genotyping** runs via `nf-core/modules/*` wrappers:
  - `trust4`: de novo TCR/BCR CDR3 assembly from bulk or single-cell RNA-seq BAM/FASTQ (`*_report.tsv`, `*_airr.tsv`).
  - `optitype`: ILP-based 4-digit MHC Class I typing from DNA/RNA BAMs (`*_result.tsv`).
  - `arcashla`: RNA-seq BAM HLA Class I/II genotyping (`*.genotype.json`).
  - `seq2hla`: RNA-seq FASTQ Class I/II HLA typing and locus expression (`*HLAgenotype4digits`).
  - `hlala`: graph-based Class I/II HLA typing from WGS/WES/long reads (`R1_bestguess_G.txt`).
- **Downstream statistical & translational analysis** runs deterministically in `plugin:imm` (`imm_*`).

## 1. Repertoire Harmonization & QC (`imm_repertoire_inspect`)

- **Supported Formats**:
  - `airr`: AIRR Rearrangement TSV/CSV (`junction_aa`/`cdr3_aa`, `v_call`, `j_call`, `locus`, `productive`, `duplicate_count`/`consensus_count`, `cell_id`, `sample_id`).
  - `cellranger_vdj`: 10x Genomics V(D)J `filtered_contig_annotations.csv` (`barcode`, `chain`, `v_gene`, `j_gene`, `cdr3`, `cdr3_nt`, `productive`, `umis`, `reads`).
  - `trust4`: TRUST4 `*_report.tsv` (`#count`/`count`, `frequency`, `CDR3nt`, `CDR3aa`, `V`, `D`, `J`, `C`) or `*_airr.tsv`.
  - `mixcr`: MiXCR `clones.tsv` (`cloneCount`/`readCount`, `cloneFraction`, `aaSeqCDR3`, `nSeqCDR3`, `allVHitsWithScore`/`bestVHit`, `allJHitsWithScore`/`bestJHit`).
- **QC Metrics**:
  - Productive %, out-of-frame %, and stop-codon (`*`) % computed prior to optional `--productive_only` filtering.
  - Chain composition (`TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGK`, `IGL`).
  - Mean/SD CDR3 AA length and **Kyte-Doolittle (1982) hydropathy (GRAVY)** across productive CDR3 sequences.

## 2. Diversity & Clonal Architecture (`imm_diversity`)

For a sample with $S$ unique productive clonotypes of counts $c_1 \ge c_2 \ge \dots \ge c_S > 0$, total $N = \sum_{i=1}^S c_i$, and relative frequencies $p_i = c_i / N$:

- **Shannon Entropy**:
  $$H_{\text{nats}} = -\sum_{i=1}^S p_i \ln p_i, \qquad H_{\text{bits}} = \frac{H_{\text{nats}}}{\ln 2}$$
- **Pielou Clonality**:
  $$\text{Clonality} = 1 - \frac{H_{\text{nats}}}{\ln S} \in [0, 1] \quad (S > 1)$$
  Higher clonality denotes stronger oligoclonal expansion; $0$ denotes uniform evenness.
- **Simpson Indices**:
  - Gini-Simpson: $1 - \sum_{i=1}^S p_i^2$
  - Inverse Simpson: $\left(\sum_{i=1}^S p_i^2\right)^{-1}$
- **$D_{50}$**: minimum fraction (%) of top-ranked clonotypes accounting for $\ge 50\%$ of total reads/UMIs:
  $$D_{50} = \frac{100 \cdot k^*}{S}, \quad k^* = \min\left\{k : \sum_{i=1}^k p_i \ge 0.50\right\}$$
- **Clonal Homeostasis Space**:
  - Hyperexpanded: $p_i \ge 1\%$
  - Large: $0.1\% \le p_i < 1\%$
  - Medium: $0.01\% \le p_i < 0.1\%$
  - Small: $p_i < 0.01\%$
- **Group Statistics**: Two-group two-sided Mann-Whitney U test or $\ge 3$-group Kruskal-Wallis H-test with Benjamini-Hochberg FDR across tested diversity metrics.

## 3. Clonotype Overlap & Longitudinal Tracking (`imm_clonotype_tracking`)

Clonotypes are keyed by `CDR3AA|V_gene` (falling back to `CDR3AA` if V gene is unavailable).

- **Jaccard Index**:
  $$J(A, B) = \frac{|K_A \cap K_B|}{|K_A \cup K_B|}$$
- **Morisita-Horn Similarity**:
  $$C_{MH}(A, B) = \frac{2 \sum_k p_{k,A} p_{k,B}}{\sum_k p_{k,A}^2 + \sum_k p_{k,B}^2}$$
- **Public Clonotypes & 4-mer Motifs**: Identifies clonotypes shared across $\ge 2$ samples and evaluates contiguous 4-mer CDR3 AA motif enrichment in shared/expanded clonotypes vs singleton clonotypes via Fisher's exact test + BH FDR.
- **Longitudinal Expansion/Contraction**: Between `baseline_sample` ($N_0$) and `followup_sample` ($N_1$), each clonotype is tested via two-sided Fisher's exact test on $\begin{pmatrix} c_{k,1} & N_1 - c_{k,1} \\ c_{k,0} & N_0 - c_{k,0} \end{pmatrix}$ with BH FDR (`expanded` if $\text{FDR} < 0.05$ and $p_{k,1} > p_{k,0}$; `contracted` if $\text{FDR} < 0.05$ and $p_{k,1} < p_{k,0}$).

## 4. HLA Consensus, Zygosity & Approximate HED (`imm_hla_summary`)

- **2-Field Standardization**: Normalizes calls from OptiType, arcasHLA, seq2HLA, HLA-LA, and clinical cohort tables to `HLA-X*XX:XX` across Class I (`HLA-A`, `HLA-B`, `HLA-C`) and Class II (`HLA-DRB1`, `HLA-DQB1`, `HLA-DPB1`). Multi-caller inputs for the same sample are reconciled by consensus voting.
- **Zygosity & Immune Escape**: Homozygosity at Class I loci restricts peptide presentation diversity and is flagged (`mhc_i_homozygosity_flag`).
- **Approximate HED (Pocket-Grantham)**: Computes pairwise Grantham (1974) physicochemical distance across curated B/C/D/E/F peptide-binding groove residue profiles (plus a calibrated 2nd-field micro-divergence term for distinct alleles within the same allele group). Identical 2-field alleles strictly yield $\text{HED} = 0.0$.

## 5. Neoantigen Prioritization (`imm_neoantigen_prioritize`)

- **Binding Affinity**:
  - Uses upstream `ic50_mt` / `ic50_wt` directly when present (`strong_binder <= 50 nM`, `weak_binder <= ic50_cutoff_nm`, default `500 nM`).
  - When `ic50_mt` is absent, estimates 9-mer affinity across overlapping mutant windows via a curated HLA Class I anchor-motif PSSM (`P2` and `C-terminus` primary anchors + secondary contacts) and explicitly flags `ic50_source = "pssm_estimate"`.
- **Differential Agretopicity Index (DAI)**:
  $$\text{DAI} = \frac{\text{IC}_{50,\text{WT}}}{\max(\text{IC}_{50,\text{MT}}, 10^{-3})}$$
- **Anchor vs. TCR-Contact Mutation**: Classifies equal-length epitopes by whether substitutions hit binding anchors (`P2` or `C-term`) vs TCR-facing residues (`P3..P_{L-2}`).
- **Expression, Clonality & Composite Immunogenicity Score (`0..100`)**: Applies `ic50_mt_nm <= ic50_cutoff_nm`, `tpm >= min_tpm`, and clonality (`ccf >= 0.5` or `vaf >= 0.1`) gates, ranking candidates by `immunogenicity_score` combining log-IC50 (40 pts), DAI/anchor-vs-TCR balance (25 pts), foreignness & hydropathy (15 pts), and expression/clonality (20 pts).
