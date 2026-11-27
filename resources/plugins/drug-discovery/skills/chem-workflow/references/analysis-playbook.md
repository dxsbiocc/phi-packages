# Small-Molecule Drug Discovery & Cheminformatics Playbook

## Dual-Engine Architecture (`rdkit` vs `pure_python_smiles`)

Every script in `chem-workflow` reports `"backend": "rdkit"` when RDKit is importable (`plugin:chem` runtime environment) and automatically falls back to `"backend": "pure_python_smiles"` in minimal Python environments lacking RDKit:

- **RDKit Engine**: `MolStandardize` (`LargestFragmentChooser`, `Uncharger`), canonical isomeric SMILES, `InChIKey`, `Crippen.MolLogP`, `rdMolDescriptors.CalcTPSA`, `QED.qed`, `FilterCatalog` (`PAINS` + `BRENK`), `MurckoScaffold`, and 1024-bit radius-2 Morgan fingerprints (`rdFingerprintGenerator` / `AllChem`).
- **Pure-Python SMILES Graph Fallback**: Parses organic SMILES (`B,C,N,O,P,S,F,Cl,Br,I,c,n,o,s,p`), bracket atoms (`[...]`), chiral/E-Z annotations, branches, and `%nn` ring closures; strips salts (retaining largest heavy-atom component); neutralizes simple protonated amines/carboxylates/phenolates; computes canonical symmetry-refined graph SMILES + deterministic `PYKEY-...` structure keys; evaluates Wildman-Crippen-style `cLogP`, Ertl `TPSA`, Bickerton `QED`, ring/scaffold pruning, 1024-bit circular radius-2 atom-neighborhood fingerprints, and 15+ classic PAINS/Brenk toxicophores.

---

## 1. Standardization (`chem_standardize`)

- **Salt/Solvent Stripping**: Split disconnected `.` components and retain the component with the maximum heavy-atom count (tie-breaking on carbon count).
- **Neutralization**: Deprotonate simple conjugate acids (`[NH3+]`, `[NH2+]`, `[NH+]`) and protonate conjugate bases (`[O-]C=O`, phenolate `[O-]c`) while preserving quaternary ammonium centers (`[N+](C)(C)(C)C`) and zwitterionic nitro groups (`[N+](=O)[O-]`).
- **Deduplication**: Collapse identical canonical structure keys (`InChIKey` or `PYKEY`) unless `--keep_duplicates` is passed.

## 2. Medicinal Chemistry Triage (`chem_medchem_triage`)

- **Lipinski Rule of Five (Ro5)**: Counts violations among `MW > 500`, `cLogP > 5`, `HBD > 5`, `HBA > 10`.
- **Veber Oral Bioavailability**: Passes when `RotBonds <= 10` and `TPSA <= 140 Å²`.
- **QED (Bickerton et al., 2012)**: Weighted geometric mean of desirability functions across `MW`, `cLogP`, `HBA`, `HBD`, `TPSA`, `RotBonds`, `AromaticRings`, and `Alerts` in `[0, 1]`.
- **Synthetic Complexity (`1..10`)**: Ring/spiro/bridged/stereocenter/heteroatom complexity heuristic normalized to `[1.0, 10.0]` (lower = easier synthesis).
- **Taylor-Butina Clustering**: Greedy neighbor-count leader clustering on 1024-bit Morgan/circular Tanimoto distance `d = 1 - sim <= cluster_cutoff` (default `0.4`, i.e., `sim >= 0.6`).
- **Triage Tiers**:
  - `alert_flagged`: $\ge 1$ PAINS or Brenk alert
  - `beyond_ro5`: no alerts and $\ge 2$ Lipinski violations
  - `lead_like`: no alerts, `MW <= 350`, `cLogP <= 3.5`, `RotBonds <= 7`, `QED >= 0.5`
  - `drug_like`: remaining clean compounds (`<= 1` Ro5 violation)

## 3. SAR & Activity Cliffs (`chem_sar_cliff`)

- **Potency Normalization**:
  - `nM`: $\text{pIC}_{50} = 9 - \log_{10}(\text{IC}_{50,\text{nM}})$
  - `uM`: $\text{pIC}_{50} = 6 - \log_{10}(\text{IC}_{50,\mu\text{M}})$
- **Efficiency Metrics**:
  - **Ligand Efficiency**: $\text{LE} = \frac{1.37 \cdot \text{pIC}_{50}}{N_{\text{heavy}}}$ (kcal/mol per heavy atom)
  - **Lipophilic Ligand Efficiency**: $\text{LLE} = \text{pIC}_{50} - \text{cLogP}$
- **Guha–Van Drie SALI & Activity Cliffs**:
  $$\text{SALI}_{ij} = \frac{|\text{pIC}_{50,i} - \text{pIC}_{50,j}|}{1 - \text{sim}_{ij} + 0.001}$$
  Pairs with $\text{sim}_{ij} \ge \text{similarity\_threshold}$ (default `0.55`) and $|\Delta\text{pIC}_{50}| \ge 1.0$ ($\ge 10\times$ potency jump) are flagged as **Activity Cliffs**.

## 4. Physicochemical ADMET Profiling (`chem_admet_profile`)

- **BOILED-Egg (Daina & Zoete, 2016)**: Elliptical classification in `(TPSA, cLogP)` space distinguishing `BBB+` (yolk), `HIA+ (BBB-)` (white), and `Poorly Absorbed` (outside).
- **Wager CNS MPO (0–6)**: Composite desirability score summing `cLogP`, `cLogD`, `MW`, `TPSA`, `HBD`, and basic-amine `pKa` proxy (`>= 4.0` favorable for CNS penetration).
- **Delaney ESOL**: $\log S = 0.16 - 0.63\,\text{cLogP} - 0.0062\,\text{MW} + 0.066\,\text{RotBonds} - 0.74\,\text{AP}$ (`AP` = aromatic heavy-atom fraction).
- **Safety Heuristics**:
  - **Pfizer 3/75**: `cLogP > 3.0` and `TPSA < 75 Å²`
  - **GSK 4/400**: `MW > 400` and `cLogP > 4.0`
  - **hERG Liability**: `has_basic_amine` and `cLogP > 3.7`

## 5. Docking & Co-Folding Evaluation (`chem_docking_eval`)

- Converts raw score to favorable binding magnitude `fav_score` (`-score` when `higher_is_better=false`, `+score` when `higher_is_better=true`).
- **Size Normalization**:
  - $\text{LE}_{\text{dock}} = \text{fav\_score} / N_{\text{heavy}}$
  - **SILE (Nissink, 2009)**: $\text{SILE} = \text{fav\_score} / N_{\text{heavy}}^{0.3}$
- Filters PAINS/Brenk alerts before selecting `top_n` cluster-diverse representatives (`Tanimoto < 0.65` diversity pass with backfill if needed), supporting optional `aux_score_col` consensus Z-score and rank fusion.
