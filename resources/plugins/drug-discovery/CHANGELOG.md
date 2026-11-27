# Changelog — drug-discovery

## 0.1.0 (2026-11-27)

Initial release of the small-molecule drug discovery, ADMET, and virtual-screening plugin (`toolPrefix: chem`, agent `DrugDiscovery`, private locked environment `plugin:chem`).

- `chem_standardize`: Salt/solvent stripping, functional-group neutralization, canonical SMILES, `InChIKey`/graph-hash deduplication, and stereochemistry/validity QC.
- `chem_medchem_triage`: QED, synthetic complexity, Lipinski Ro5 & Veber rules, PAINS/Brenk alerts, Bemis-Murcko scaffolds, and Taylor-Butina ECFP4 clustering.
- `chem_sar_cliff`: `pIC50`, `LE`, `LLE`, pairwise ECFP4 Tanimoto similarity, Guha-Van Drie `SALI`, and activity-cliff detection.
- `chem_admet_profile`: BOILED-Egg (`HIA`/`BBB`), Wager CNS MPO, Delaney ESOL `logS`, Pfizer `3/75`, GSK `4/400`, and hERG liability profiling.
- `chem_docking_eval`: PAINS filtering, size-normalized `LE_dock` & `SILE`, optional consensus scoring, and cluster-diverse virtual-screening hit selection.
