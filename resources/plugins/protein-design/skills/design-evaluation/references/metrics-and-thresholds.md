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

## 4. Pre-Wetlab Developability Audit (`pd_developability_audit`)

| Metric / Liability | Target / Safe Range | Warning / Risk Rule |
| --- | --- | --- |
| Isoelectric point (`pI`) | `|pI - buffer_ph| >= 1.0` | Flag `pi_buffer_clash` if `|pI - buffer_ph| < 1.0` (isoelectric precipitation risk) |
| UV280 extinction ($\varepsilon_{280}$, $A_{280}^{0.1\%}$) | `n_W + n_Y >= 1` | Flag `uv280_blind` if $n_W = 0$ and $n_Y = 0$ (NanoDrop A280 blind; add Trp or use BCA) |
| Guruprasad Instability Index | `<= 40.0` | `> 40.0` flagged unstable (`> 45.0` escalates risk) |
| 9-mer hydrophobic patch (KD) | `<= 2.2` | `> max_hydrophobic_patch` (default `2.2`) flags surface aggregation risk |
| Unpaired Cysteines | `n_cys % 2 == 0` | Odd `n_cys` flags `unpaired_cys` (thiol scrambling / covalent dimerization) |
| Chemical hotspots | `NG/NS/NT`, `DG/DS`, `N-X-[ST]` | Counted and localized for deamidation, isomerization, and mammalian N-glycosylation |
| Overall tier | `pass` | `high_risk` (severe flags: `uv280_blind`, `unpaired_cys`, `pi_buffer_clash`, `hydrophobic_patch_risk`, or `II > 45`), `caution` (moderate warnings), `pass` (clean) |
