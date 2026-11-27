# Molecular Cloning Analysis Playbook

## 1. Primer Thermodynamics (`mc_design_primers`)

### SantaLucia (1998) Nearest-Neighbor + Owczarzy et al. (2008) Salt Correction
- **Unified NN $\Delta H^\circ, \Delta S^\circ$**: Sum dinucleotide steps (`AA/TT`, `AT/TA`, `TA/AT`, `CA/GT`, `GT/CA`, `CT/GA`, `GA/CT`, `CG/GC`, `GC/CG`, `GG/CC`) plus terminal initiation (`G/C`: $+0.1\text{ kcal/mol}, -2.8\text{ cal/(mol}\cdot\text{K)}$; `A/T`: $+2.3\text{ kcal/mol}, +4.1\text{ cal/(mol}\cdot\text{K)}$).
- **1 M $\text{Na}^+$ baseline**: $T_m(1\text{M}) = \frac{1000 \cdot \Delta H^\circ}{\Delta S^\circ + R \ln(C_T / 4)}$.
- **Owczarzy (2008) divalent $\text{Mg}^{2+}$ / dNTP / $\text{Na}^+$ regimes**:
  - Free $\text{Mg}^{2+}$: $[\text{Mg}^{2+}]_{\text{free}} = \max(0, [\text{Mg}^{2+}] - [\text{dNTP}])$.
  - Using $R = \sqrt{[\text{Mg}^{2+}]_{\text{free}}} / [\text{Na}^+]$:
    - $R < 0.22$: monovalent-dominated correction.
    - $0.22 \le R < 6.0$: mixed $\text{Na}^+$/$\text{Mg}^{2+}$ correction (Eq. 16–20).
    - $R \ge 6.0$: divalent $\text{Mg}^{2+}$-dominated correction.
- **Secondary structure & mispriming**: Score GC% (`40–60%`), 3' GC clamp (`1–2` G/C in last 5 bases; penalize `>=4`), homopolymers (`>=4`), hairpin stem stability, self-dimer & 3'-anchored cross-dimer $\Delta G_{37}$, $|\Delta T_m| \le 2.5^\circ\text{C}$, and optional 3'-12mer hits in `background_fasta`.

## 2. Host Codon Optimization (`mc_codon_optimize`)

- Supports `E_coli_K12`, `HEK293_Human`, `CHO_Cricetulus`, and `S_cerevisiae`.
- Maximizes Codon Adaptation Index ($\text{CAI} = (\prod w_i)^{1/L}$) while enforcing:
  1. **100% identical translated amino-acid sequence** (hard assertion).
  2. **0 forbidden restriction sites** (`BsaI`, `BsmBI`, `BbsI`, `EcoRI`, `BamHI`, `HindIII`, `NdeI`, `XhoI`, or user-specified) on both strands.
  3. **0 homopolymers** (`TTTT`, `AAAAA`, `GGGGGG`, `CCCCCC`) and **0 cryptic Shine-Dalgarno** (`AGGAGG`) in *E. coli*.
  4. **30-bp sliding-window GC** inside `[gc_window_min, gc_window_max]` and relaxed 5' GC across codons 2–8.

## 3. Construct Assembly (`mc_assemble_construct`)

- **GoldenGate**: Digests linear parts or circular vectors with `BsaI` (`GGTCTC`), `BsmBI`/`Esp3I` (`CGTCTC`), `BbsI`/`BpiI` (`GAAGAC`), or `SapI` (`GCTCTTC`), builds the directed overhang graph, and flags palindromic overhangs (`oh == revcomp(oh)`), duplicate/revcomp collisions, 1-bp mismatch crosstalk, and retained internal Type IIS sites.
- **Gibson**: Detects 15–60 bp terminal homologous overlaps (including circular wrap-around), computes overlap $T_m$, and warns on short (`<15 bp`) or low-$T_m$ (`<48 °C`) junctions.
- **RestrictionLigation**: Supports 35+ Type II enzymes including compatible cohesive ends (`BamHI`/`BglII`, `XbaI`/`SpeI`/`NheI`/`AvrII`, `SalI`/`XhoI`).

## 4. Sanger ABIF `.ab1` Verification (`mc_sanger_verify`)

- Parses binary `ABIF` directory tags (`FWO_`, `PBAS`, `PCON`, `PLOC`, `DATA9..12`).
- Applies **Mott quality trimming** ($p_i = 10^{-Q_i/10}$, cutoff `trim_error_prob = 0.05`).
- Auto-orients strands and supports **circular plasmid origin wrap-around** alignment.
- Reports SNVs, insertions, deletions, frameshift flags, and `.ab1` secondary peak ratios $I_2/(I_1+I_2) \ge 0.15$.
