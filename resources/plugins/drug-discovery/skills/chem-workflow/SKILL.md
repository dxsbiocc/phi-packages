---
name: chem-workflow
description: >-
  Small-molecule drug discovery and cheminformatics workflow: SMILES salt
  stripping, neutralization, canonicalization, and deduplication; Lipinski Ro5,
  Veber, QED, synthetic complexity, PAINS/Brenk alerts, Bemis-Murcko scaffolds,
  and Taylor-Butina clustering; SAR potency normalization (pIC50), LE/LLE, and
  SALI activity-cliff mapping; BOILED-Egg, CNS MPO, Delaney ESOL, Pfizer 3/75,
  GSK 4/400, and hERG ADMET profiling; and size-normalized docking/co-folding
  hit evaluation (LE_dock, SILE, consensus ranking, diversity selection). Use
  for SMILES, RDKit, medchem triage, QED, PAINS, Butina, SAR, activity cliff,
  SALI, ADMET, BOILED-Egg, ESOL, CNS MPO, hERG, docking, or virtual screening.
phi:
  environment: plugin:chem
  attachTo:
    - DrugDiscovery
  scripts:
    - name: standardize
      description: >-
        Standardize small-molecule structures (TSV/CSV/SMI/SDF): strip salts and
        solvents keeping the largest organic fragment, neutralize ionized groups,
        canonicalize SMILES, compute structure keys (InChIKey or graph hash),
        flag invalid SMILES and stereochemistry, and deduplicate compounds.
      run: [python, ./scripts/chem_standardize.py]
      args:
        type: object
        additionalProperties: false
        required: [compounds, output_dir]
        properties:
          compounds:
            type: string
            format: input-path
            description: Input TSV/CSV/SMI/SDF compound file containing SMILES and optional compound IDs.
          output_dir:
            type: string
            format: project-path
            description: Output directory for standardized_compounds.tsv and standardization_summary.png.
          smiles_col:
            type: string
            description: Optional column name containing SMILES strings (auto-detected if omitted).
          id_col:
            type: string
            description: Optional column name containing compound identifiers (auto-detected if omitted).
          keep_duplicates:
            type: boolean
            description: Retain duplicate canonical structures instead of collapsing to unique keys (default false).
      approval: write
      output: ./schemas/standardize-result.json
      timeoutSeconds: 600
    - name: medchem_triage
      description: >-
        Medicinal chemistry library triage: compute MW, cLogP, TPSA, HBD, HBA,
        RotBonds, Fsp3, QED, synthetic complexity (1-10), Lipinski Ro5 & Veber
        compliance, PAINS & Brenk structural alerts, Bemis-Murcko scaffolds, and
        Taylor-Butina ECFP4 diversity clusters with lead-like/drug-like tiers.
      run: [python, ./scripts/chem_medchem_triage.py]
      args:
        type: object
        additionalProperties: false
        required: [compounds, output_dir]
        properties:
          compounds:
            type: string
            format: input-path
            description: Input TSV/CSV table containing compound SMILES and optional IDs.
          output_dir:
            type: string
            format: project-path
            description: Output directory for medchem_triage.tsv, scaffold_summary.tsv, and chemical_space_qed.png.
          smiles_col:
            type: string
            description: Optional SMILES column name (auto-detected if omitted).
          id_col:
            type: string
            description: Optional compound ID column name (auto-detected if omitted).
          cluster_cutoff:
            type: number
            description: Taylor-Butina Tanimoto distance cutoff in (0, 1) (default 0.4, corresponding to similarity >= 0.6).
      approval: write
      output: ./schemas/medchem-triage-result.json
      timeoutSeconds: 600
    - name: sar_cliff
      description: >-
        Structure-Activity Relationship (SAR) and activity-cliff analysis:
        normalize potency to pIC50, compute Ligand Efficiency (LE) and
        Lipophilic Ligand Efficiency (LLE), build pairwise ECFP4 Tanimoto
        similarity matrix, compute Guha-Van Drie SALI, and flag Activity Cliffs.
      run: [python, ./scripts/chem_sar_cliff.py]
      args:
        type: object
        additionalProperties: false
        required: [compounds, potency_col, output_dir]
        properties:
          compounds:
            type: string
            format: input-path
            description: Input TSV/CSV compound table with SMILES and measured potency values.
          potency_col:
            type: string
            description: Column name containing potency measurements (e.g. ic50_nm, ki_nm, or pic50).
          output_dir:
            type: string
            format: project-path
            description: Output directory for sar_metrics.tsv, activity_cliffs.tsv, and sar_lle_cliff_landscape.png.
          potency_scale:
            type: string
            enum: [nM, uM, pIC50]
            description: Scale of potency_col values (default nM).
          smiles_col:
            type: string
            description: Optional SMILES column name (auto-detected if omitted).
          id_col:
            type: string
            description: Optional compound ID column name (auto-detected if omitted).
          similarity_threshold:
            type: number
            description: Minimum ECFP4 Tanimoto similarity to flag an activity cliff pair (default 0.55).
      approval: write
      output: ./schemas/sar-cliff-result.json
      timeoutSeconds: 600
    - name: admet_profile
      description: >-
        Physicochemical ADMET and safety liability profiling: BOILED-Egg HIA and
        BBB permeation classification, Wager CNS MPO (0-6), Delaney ESOL aqueous
        solubility logS, Pfizer 3/75 promiscuity risk, GSK 4/400 ADMET risk, and
        hERG basic-amine lipophilicity liability flagging.
      run: [python, ./scripts/chem_admet_profile.py]
      args:
        type: object
        additionalProperties: false
        required: [compounds, output_dir]
        properties:
          compounds:
            type: string
            format: input-path
            description: Input TSV/CSV compound table with SMILES and optional compound IDs.
          output_dir:
            type: string
            format: project-path
            description: Output directory for admet_profile.tsv and boiled_egg_admet.png.
          smiles_col:
            type: string
            description: Optional SMILES column name (auto-detected if omitted).
          id_col:
            type: string
            description: Optional compound ID column name (auto-detected if omitted).
      approval: write
      output: ./schemas/admet-profile-result.json
      timeoutSeconds: 600
    - name: docking_eval
      description: >-
        Virtual-screening docking and co-folding score evaluation: join ligand
        structures, filter PAINS/Brenk alerts, compute size-normalized LE_dock
        and SILE, fuse optional auxiliary scores via consensus ranking, and
        select cluster-diverse top hits.
      run: [python, ./scripts/chem_docking_eval.py]
      args:
        type: object
        additionalProperties: false
        required: [docking_scores, score_col, output_dir]
        properties:
          docking_scores:
            type: string
            format: input-path
            description: Docking/affinity score TSV/CSV (Vina, GNINA, Glide, or boltz2) with compound ID, score, and optional SMILES.
          score_col:
            type: string
            description: Primary docking or affinity score column name.
          output_dir:
            type: string
            format: project-path
            description: Output directory for docking_hit_list.tsv and docking_le_pareto.png.
          compounds:
            type: string
            format: input-path
            description: Optional TSV/CSV table mapping compound IDs to SMILES strings.
          higher_is_better:
            type: boolean
            description: Set true when larger primary scores are better (e.g. CNNscore or pKd); default false for negative deltaG/Vina kcal/mol.
          aux_score_col:
            type: string
            description: Optional secondary score column (e.g. cnn_score) for consensus Z-score and rank fusion.
          top_n:
            type: integer
            description: Maximum number of cluster-diverse hits to select (default 25).
      approval: write
      output: ./schemas/docking-eval-result.json
      timeoutSeconds: 600
---

# Small-Molecule Drug Discovery, ADMET & Virtual Screening Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for the mathematical formulations, RDKit vs pure-Python SMILES graph engine behavior, ADMET alert rules, and docking efficiency formulas.

## Tool Overview

```
Raw SMILES / SDF / SMI     ──chem_standardize────▶ Salt-stripped, neutralized, deduplicated compounds
Standardized library       ──chem_medchem_triage─▶ QED, Ro5/Veber, PAINS/Brenk, Murcko, Butina clusters
Compounds + IC50/Ki/pIC50  ──chem_sar_cliff──────▶ pIC50, LE, LLE, Tanimoto matrix, SALI activity cliffs
Compounds                  ──chem_admet_profile──▶ BOILED-Egg (HIA/BBB), CNS MPO, ESOL, 3/75, 4/400, hERG
Docking / boltz2 scores    ──chem_docking_eval───▶ PAINS filter, LE_dock & SILE, consensus & diverse hits
```

## Hard Red Lines

1. **Always standardize before computing descriptors**: salt/solvent fragments and protonation states distort `MW`, `cLogP`, `TPSA`, and duplicate counts.
2. **Never silently ignore PAINS/Brenk alerts**: reactive and pan-assay interference chemotypes must remain explicitly flagged in both medchem triage and docking hit lists.
3. **Never rank docking hits solely by raw binding energy**: always inspect size-normalized `LE_dock` and `SILE` alongside scaffold diversity to avoid molecular-weight inflation bias.
4. **Treat ADMET rules as early filters**: BOILED-Egg, ESOL, CNS MPO, Pfizer `3/75`, GSK `4/400`, and hERG basic-amine alerts are prioritization heuristics, not experimental DMPK or toxicology measurements.
