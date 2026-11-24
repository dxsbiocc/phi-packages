# Metrics & Filtering Thresholds Reference

## 1. Core Starting Thresholds (PPI & Small-Molecule Binders)

| Metric | Column name(s) | Typical passing starting point | Notes |
| --- | --- | --- | --- |
| Backbone self-consistency | `ca_rmsd`, `backbone_rmsd` | `<= 2.0` Å (`<= 1.5` Å stringent) | Lower is better; measures agreement between designed backbone and refolded structure |
| Structure confidence | `complex_plddt` (0–1) / `plddt` (0–100) | `complex_plddt >= 0.80` (`plddt >= 80`) | Both Boltz-2 and AlphaFold 3 emit normalized `complex_plddt` in `0..1` |
| Interface confidence | `iptm` (or `ligand_iptm`) | `>= 0.70` (`>= 0.80` stringent) | Primary predictor of complex interface plausibility |
| Interface PAE (AF3) | `min_pae_interaction` | `<= 5.0` Å (`<= 10.0` Å permissive) | Off-diagonal minimum of `chain_pair_pae_min` |
| Steric quality | `clash_count` | `== 0` | Non-bonded inter-residue heavy-atom pairs < 2.0 Å |
| Interface packing | `interface_contacts`, `interface_hbonds` | `interface_contacts >= 15`, `interface_hbonds >= 2` | Verify binder actually contacts target chain rather than folding in isolation |
| Boltz-2 ligand hit probability | `affinity_probability_binary` | `>= 0.5` (hit discovery) | Boltz-2 only; separates binders from decoys |
| Boltz-2 ligand affinity | `affinity_pred_value` (`log10(IC50 uM)`) | Lower is stronger (`-3` ≈ 1 nM, `0` ≈ 1 uM) | Use to rank active binders, not to replace binary hit filtering |

## 2. Example `DESIGN_FILTER` Boolean DSL Expressions

- **Boltz-2 PPI binder**:
  `complex_plddt >= 0.80 and iptm >= 0.70 and ca_rmsd <= 2.0 and interface_contacts >= 15 and clash_count == 0`
- **AlphaFold 3 PPI binder**:
  `complex_plddt >= 0.80 and iptm >= 0.75 and min_pae_interaction <= 5.0 and ca_rmsd <= 2.0 and clash_count == 0`
- **Boltz-2 small-molecule binder**:
  `complex_plddt >= 0.80 and ligand_iptm >= 0.75 and affinity_probability_binary >= 0.5 and ligand_contacts >= 10 and clash_count == 0`

## 3. Red Lines

- Never rank candidates pooled from different structure predictors (`boltz2` vs `alphafold3`) on raw `confidence_score` without separating or re-folding all top hits with one common oracle.
- Always inspect `unique_parents_selected` in `pd_funnel_report` — ordering 24 sequences that all derive from a single backbone carries high experimental correlation risk.
