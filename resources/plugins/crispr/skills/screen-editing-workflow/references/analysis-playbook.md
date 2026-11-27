# Pooled Screen & Post-Editing Analysis Playbook

## 1. Upstream Wrappers vs. Local Tools

- **`nf-core/modules/mageck-count`**: converts raw pooled screen FASTQ reads into `.count.txt` sgRNA count matrices.
- **`nf-core/modules/crisprcleanr-normalize`**: removes copy-number-amplification false positives in cancer cell-line viability screens prior to hit calling.
- **`nf-core/modules/mageck-test` / `mageck-mle`**: optional upstream MAGeCK CLI alternatives for multi-arm MLE design matrices.
- **`nf-core/modules/crispresso2`**: demultiplexes and aligns raw deep-sequencing amplicon FASTQ reads to produce `Alleles_frequency_table_around_sgRNA_*.txt`.

## 2. Pooled Screen QC (`crispr_screen_qc`)

- **Gini index (Lorenz inequality)**:
  $$G = \frac{\sum_{i=1}^n (2i - n - 1) x_{(i)}}{n \sum_{i=1}^n x_{(i)}}$$
  Baseline/plasmid libraries should have $G \le 0.15$; $G > 0.25$ indicates representation skew/bottlenecking and $G > 0.30$ flags severe library collapse.
- **NNMD (Null-Normalized Mean Difference)** against a reference/T0 baseline:
  $$\text{NNMD} = \frac{\mu_{\text{essential}} - \mu_{\text{nonessential}}}{\sigma_{\text{nonessential}}}$$
  Strong essentiality separation achieves $\text{NNMD} < -1.25$ (and essentiality ROC-AUC $> 0.85$).

## 3. Gene-Level Hit Calling (`crispr_screen_hits`)

1. **Count normalization**: DESeq median-of-ratios (`median`), non-targeting control median-of-ratios (`control_sgrna`), or total read scaling (`total`).
2. **Moderated sgRNA variance & $z$-scores**: estimates empirical mean-variance trend across sgRNAs on normalized counts and shrinks per-sgRNA variance toward the trend + Poisson floor, yielding signed sgRNA $z$-scores and $\log_2\text{FC}$.
3. **$\alpha$-RRA (Li et al., 2014)**: ranks all $M$ sgRNAs separately for depletion and enrichment ($u_i = \text{rank}_i / M$). For each gene with $k$ sgRNAs and sorted uniform ranks $u_{(1)} \le \dots \le u_{(k)}$:
   $$\rho = \min_{1 \le j \le k:\, u_{(j)} \le \alpha} I_{u_{(j)}}(j,\, k - j + 1)$$
   Exact gene-level $p$-values are calibrated against `n_permutations` (default `2000`) seeded null draws stratified by sgRNA count $k$, followed by **Benjamini–Hochberg FDR** (`neg_fdr`, `pos_fdr`) and Stouffer's combined $Z$.

## 4. Post-Editing Outcome Quantification (`crispr_editing_quant`)

- Quantifies reads inside `[cut - window_bp, cut + window_bp)` (and across protospacer positions `1..L`) via banded/global Needleman-Wunsch alignment or CRISPResso2 aligned pairs:
  - **Unmodified**: matches WT across window.
  - **Frameshift Indel**: `n_ins + n_del > 0` with `|n_ins - n_del| % 3 != 0` (productive KO).
  - **In-Frame Indel**: `n_ins + n_del > 0` with `|n_ins - n_del| % 3 == 0`.
  - **Substitution Only**: `n_ins == 0, n_del == 0, n_sub > 0` (expected pure base-editing product).
  - **HDR / Imperfect HDR**: matches `hdr_seq` signature edits cleanly vs accompanied by unintended indels.
- For `CBE` (`C->T`) and `ABE` (`A->G`), reports base conversion across protospacer positions `1..L`, target conversion % inside `[4..8]`, bystander edit %, and indel byproduct %.

## 5. Chemogenetic Drug-vs-Vehicle Screening (`crispr_drugz`)

- **When to prefer `crispr_drugz` over `crispr_screen_hits`**:
  - In chemogenetic modifier screens comparing **Drug-treated vs. DMSO-treated** arms at the same endpoint (`IC20`–`IC30` sub-lethal drug dose), low-abundance sgRNAs exhibit inflated fold-change variance. **DrugZ** (Collimore & Hart, *Genome Medicine* 2019) normalizes each sample to $10^7$ reads, orders sgRNAs by descending control abundance, fits a **monotonically non-decreasing empirical-Bayes MAD variance** across sliding windows (`window_size=500`), and aggregates replicate/guide $z$-scores into centered unit-variance **`normZ`**:
    $$\text{normZ}_g = \frac{\sum_{i \in g} z_i}{\sqrt{n_{\text{obs}, g}}} \quad (\text{re-centered and MAD-standardized across genes})$$
  - **`normZ < 0` (`fdr_synth_lethal < 0.05`)**: synthetic-lethal **drug sensitizers** (knockout enhances drug killing).
  - **`normZ > 0` (`fdr_suppressor < 0.05`)**: **drug-resistance suppressors** (knockout confers resistance, e.g., drug target, importer, or pro-drug activating enzyme).
  - Pass `--remove_genes` for depleted core-essential genes if late-timepoint pan-essential dropout distorts low-count variance.
