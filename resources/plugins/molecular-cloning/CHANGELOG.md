# Changelog — molecular-cloning

## 0.1.0 (2026-11-26)

Initial release: Molecular cloning, primer design, codon optimization, construct assembly, and Sanger ABIF verification plugin (`toolPrefix: mc`, agent `MolecularCloning`, private locked environment `plugin:mc`).

- `mc_design_primers`: `PCR`, `qPCR`, `Gibson`, `GoldenGate`, and `SiteDirectedMutagenesis` primer design using exact SantaLucia (1998) nearest-neighbor thermodynamics + Owczarzy et al. (2008) divalent $\text{Mg}^{2+}$/dNTP/$\text{Na}^+$ salt correction, GC clamp / homopolymer / hairpin / self- & cross-dimer $\Delta G_{37}$ scoring, and optional background 3'-12mer mispriming scan.
- `mc_codon_optimize`: deterministic host-specific CAI optimization (`E_coli_K12`, `HEK293_Human`, `CHO_Cricetulus`, `S_cerevisiae`) with strict 100% amino-acid identity preservation, 30-bp sliding-window GC smoothing, 5' secondary-structure relaxation, and complete elimination of forbidden restriction sites, homopolymers, and cryptic bacterial Shine-Dalgarno motifs.
- `mc_assemble_construct`: virtual `GoldenGate` (overhang graph, palindromic/duplicate/1-bp-mismatch fidelity & retained-site checks), `Gibson` (15–60 bp overlap detection & $T_m$), and `RestrictionLigation` (35+ Type II enzymes with compatible cohesive-end matching) assembly with annotated GenBank export and 2-panel construct map + virtual agarose gel.
- `mc_sanger_verify`: pure-Python binary ABIF `.ab1` chromatogram parser (`FWO_`, `PBAS`, `PCON`, `PLOC`, `DATA9..12`) + FASTQ/FASTA support, Mott Q-score trimming, circular plasmid origin wrap-around alignment, SNV/indel/frameshift calling, and secondary peak deconvolution $I_2 / (I_1 + I_2)$ for mixed colony / base-editing verification.
