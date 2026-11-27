# drug-discovery · 0.1.0

Small-molecule cheminformatics, medicinal chemistry triage, SAR activity-cliff analysis, physicochemical ADMET profiling, and virtual-screening hit evaluation plugin (`toolPrefix: chem`, agent `DrugDiscovery`) running in a locked `conda-forge` environment (`plugin:chem`: Python 3.12 + RDKit + NumPy/pandas/SciPy/Matplotlib), with a deterministic pure-Python SMILES molecular graph fallback when RDKit is unavailable.

## Tools (`chem_*`)

- `chem_standardize`: Salt/solvent stripping (retaining largest organic heavy-atom fragment), common functional-group neutralization, canonical SMILES, structure-key generation (`InChIKey` or deterministic graph hash), invalid-SMILES & stereochemistry flagging, and duplicate collapsing.
- `chem_medchem_triage`: Physicochemical descriptors (`MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotBonds`, `Fsp3`, `QED`, synthetic complexity score `1..10`), Lipinski Ro5 & Veber rules, PAINS & Brenk structural alerts, Bemis-Murcko scaffolds, Taylor-Butina ECFP4 clustering, and tiered medchem classification (`lead_like`, `drug_like`, `alert_flagged`, `beyond_ro5`).
- `chem_sar_cliff`: Potency normalization (`nM`/`uM`/`pIC50` → `pIC50`), Ligand Efficiency (`LE`), Lipophilic Ligand Efficiency (`LLE`), pairwise ECFP4 Tanimoto similarity, Guha-Van Drie `SALI`, and Activity Cliff identification (`sim >= threshold` and `|ΔpIC50| >= 1.0`).
- `chem_admet_profile`: BOILED-Egg HIA/BBB permeation model, Wager CNS MPO (`0..6`), Delaney ESOL aqueous solubility (`logS`), Pfizer `3/75` promiscuity rule, GSK `4/400` ADMET composite rule, and hERG cardiac-channel liability heuristic (`basic_amine` + `cLogP > 3.7`).
- `chem_docking_eval`: Post-docking and co-folding score evaluation (`Vina`, `GNINA`, `Glide`, `boltz2` affinity tables) with ID-based ligand joins, PAINS/Brenk filtering, size-normalized Docking Ligand Efficiency (`LE_dock`) and `SILE` (`-score / N_heavy^0.3`), optional consensus rank/Z-score fusion, and cluster-diverse hit selection.

## Scope & Division of Work with Related Modules

- **`local/boltz2`, `local/chai1`, `local/protenix`, `local/alphafold3`**: GPU structural co-folding wrappers that predict 3D protein–ligand complex coordinates and raw binding confidence/affinity tables. `drug-discovery` consumes those ligand affinity/score tables via `chem_docking_eval` to strip PAINS artifacts, correct for molecular-weight bias (`LE_dock`, `SILE`), and pick structurally diverse leads.
- **`pharmacogenomics`**: Owns raw concentration–response curve fitting (`IC50`, `EC50`, `Emax`, `AUC`) across cell lines/patients and PGx variant/star-allele associations. `drug-discovery` consumes already-fitted compound potency values (`IC50`, `Ki`, `pIC50`) in `chem_sar_cliff` to map structure–activity relationships, `LE`/`LLE`, and activity cliffs.
