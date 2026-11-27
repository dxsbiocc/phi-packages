# Guide Design, Off-Target & Library QC Playbook

## 1. Nuclease Architectures (`crispr_design_guides`)

| Nuclease | PAM | PAM Orientation | Spacer Length | Cut Site |
| --- | --- | --- | --- | --- |
| `SpCas9` | `NGG` | 3' of spacer | 20 nt | `-3` bp from PAM (pos 17) |
| `SpCas9-NG` | `NG` | 3' of spacer | 20 nt | `-3` bp from PAM (pos 17) |
| `SaCas9` | `NNGRRT` | 3' of spacer | 21 nt | `-3` bp from PAM (pos 18) |
| `Cas12a` | `TTTV` | 5' of spacer | 23 nt | Staggered `+18` / `+23` from PAM |

### Efficiency & Bae et al. (2014) Out-of-Frame MMEJ Score

1. **Sequence-context efficiency score (`0..100`)**: logistic transform of PAM-proximal nucleotide preferences (e.g., `G` favored and `C` disfavored at pos 20 for Cas9; `C` favored at cleavage pos 16/17; `TT` dinucleotide penalty in 3' seed; GC optimality in `[0.40, 0.65]`; homopolymer / `TTTT` penalties).
2. **MMEJ out-of-frame score (`0..100`)**: scans `±30 bp` flanking the cut site for microhomologies of length $k \ge 2$ separated by deletion span $d$:
   $$w(k, d) = k \cdot \exp(-d / 17)$$
   $$\text{OutOfFrameScore} = 100 \times \frac{\sum_{d \bmod 3 \neq 0} w(k, d)}{\sum_{\text{all}} w(k, d)}$$
   Higher `out_of_frame_score` predicts a higher fraction of reading-frame-disrupting deletions via MMEJ, reducing in-frame rescue alleles.
3. **Base editing (`CBE` / `ABE`)**: evaluates canonical editing window `[4..8]` (5'->3' protospacer positions). Counts editable `C` (`CBE`) or `A` (`ABE`) targets, flags secondary editable bases in `[4..8]` as **bystander edits**, and checks in-frame codon translation for premature stop creation (`CAA->TAA`, `CAG->TAG`, `CGA->TGA`, or antisense `TGG->TGA/TAG/TAA`).

### Golden Gate Cloning Oligos

- **`lentiCRISPRv2_BsmBI` / `pLKO5_BsmBI` / `pX330_BbsI`**:
  - If spacer does not start with `G`, prepend `5'-G` for U6 Pol III transcription:
    - `oligo_fwd`: `5'-CACCG` + `spacer`
    - `oligo_rev`: `5'-AAAC` + `revcomp(spacer)` + `C`
  - If spacer already starts with `G`:
    - `oligo_fwd`: `5'-CACC` + `spacer`
    - `oligo_rev`: `5'-AAAC` + `revcomp(spacer)`

## 2. Off-Target Specificity Scoring (`crispr_offtarget_scan`)

### Hsu-Zhang (2013) MIT Score
For an off-target with mismatch indices $P \subseteq \{0,\dots,19\}$ (0 = 5' PAM-distal, 19 = 3' PAM-proximal):
$$s_{\text{MIT}} = \left(\prod_{p \in P} (1 - W[p])\right) \cdot \frac{1}{\frac{19 - \bar{d}}{19}\cdot 4 + 1} \cdot \frac{1}{|P|^2} \cdot f_{\text{PAM}}$$
where $\bar{d}$ is mean pairwise distance between mismatches ($19$ if $|P|=1$) and $f_{\text{PAM}} \in \{1.0\,(\text{NGG}), 0.259\,(\text{NAG}), 0.069\,(\text{NGA})\}$.

### Doench et al. (2016) CFD Score
Product of position-and-mismatch-specific `(rN:dN)` cleavage probabilities across all mismatched positions multiplied by PAM probability (`NGG=1.0`, `NAG=0.259`, `NGA=0.069`).

### Per-Guide Specificity (`0..100`)
$$\text{Specificity} = \frac{100}{1 + \sum_{i \in \text{off-targets}} s_i}$$
where exact duplicates beyond the primary on-target locus contribute $s_i = 1.0$. Specificity `< 50` or `0` seed mismatches (`pos 13..20`) triggers an off-target warning.

> **Wrapper Escalation**: For full-genome alignment (>250 MB FASTA), invoke `nf-core/modules/bowtie` (`bowtie-build` + `bowtie-align`) or `nf-core/modules/bowtie2`.

## 3. Pooled Library Audit (`crispr_library_check`)

Before ordering oligo pools or running `nf-core/modules/mageck-count`, verify:
- Zero duplicate `sgRNA` IDs and zero unintended duplicate spacer sequences.
- Uniform spacer length and strict `ACGT` alphabet.
- `>= 3` sgRNAs per gene (typically 4–10 in genome-wide libraries) and `>= 1%` non-targeting / safe-harbor controls.
- Absence of `TTTT` Pol III terminators and cloning enzyme collisions (`BsmBI`, `BbsI`, `SapI`).

## 4. Prime Editing pegRNA & PE3/PE3b Nicking Design (`crispr_design_pegrna`)

- **SpCas9 H840A nickase geometry**: nicks the **non-target strand** (`3 bp` upstream of PAM, pos 17 of 20-nt spacer), liberating a `3'-OH` flap that hybridizes to the `3'` **Primer Binding Site (PBS)** of the pegRNA so reverse transcriptase copies the **RT Template (RTT)**.
- **Thermodynamic & structural design rules (Anzalone et al. 2019; Nelson et al. 2022)**:
  - **PBS $T_m$**: optimal around $37\text{--}42^\circ\text{C}$ (`10–15 nt`, GC `35–70%`).
  - **RTT `5'-C` penalty**: RTTs starting with cytosine (`rtt_starts_with_c: true`) base-pair with `G81` of the SpCas9 sgRNA scaffold and disrupt structure; shift `rtt_len` by `±1–2 nt` when possible.
  - **PAM/seed disruption**: edits that disrupt the `[13..20]` seed or `NGG` PAM prevent re-nicking of the edited allele; otherwise pass `cds_frame_offset` to inspect `pam_silent_mutation_hint`.
  - **PE3b vs. PE3**: **PE3b** opposite-strand `ngRNA` spacers match *only* the edited sequence (or place the PAM on the edited bases), nicking the unedited strand exclusively after prime-edited flap resolution and strongly suppressing concurrent DSB indels compared to distal `+40..+90 bp` **PE3** nicks.
  - **epegRNA (`add_epgrna_motif: true`)**: appends a linker + `37-nt tevopreQ1` RNA pseudoknot to block `3'` exonuclease degradation of the RTT/PBS extension.
