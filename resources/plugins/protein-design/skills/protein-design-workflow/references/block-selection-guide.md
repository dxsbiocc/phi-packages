# Block & Tool Selection Guide

## 1. Sequence Design (`local/modules/ligandmpnn`)

| `model_type` | Best for | Recommended `temperature` |
| --- | --- | --- |
| `protein_mpnn` | Standard protein-only backbones and PPI binders | `0.1` (production) / `0.2` (exploration) |
| `ligand_mpnn` | Small-molecule binders, enzymes, metal/cofactor complexes, nucleic-acid interfaces | `0.1`–`0.2` |
| `soluble_mpnn` | Improving solubility and expression of hydrophobic de novo scaffolds | `0.1`–`0.2` |

- Keep `output_mode: sequences` by default: downstream `boltz2` and `alphafold3` modules read `sequences.tsv` directly without needing per-sequence packed PDBs.
- Use `output_mode: structures` only when you explicitly need side-chain-packed PDBs for physics scoring or visual inspection before refolding.

## 2. Structure Prediction (`boltz2` / `alphafold3` / `protenix` / `chai1`)

| Predictor (`predictor` value) | License & weights | Highlights & key columns |
| --- | --- | --- |
| `boltz2` (`local/modules/boltz2`) | MIT code & open weights (`boltz_cache_dir`) | Fast open-license default + **ligand affinity prediction** (`affinity_pred_value`, `affinity_probability_binary`), `confidence_score`, `complex_plddt`, `iptm` |
| `alphafold3` (`local/modules/alphafold3`) | Apache-2.0 code; **non-commercial** weights (`af3_model_dir`) | Orthogonal gold-standard refolding (`ranking_score`, `complex_plddt`, `iptm`, `min_pae_interaction`, `has_clash`) |
| `protenix` (`local/modules/protenix`) | Apache-2.0 code & open weights (`protenix_cache_dir`) | Open AF3-family predictor supporting `protenix_base_default_v1.0.0` and `protenix-v2` (`ranking_score`, `complex_plddt`, `iptm`, `gpde`) |
| `chai1` (`local/modules/chai1`) | Apache-2.0 code & open weights (`chai1_cache_dir`) | Multi-modal ESM-backed folding without MSA server (`aggregate_score`, `ptm`, `iptm`, `has_inter_chain_clashes`); uses SMILES for ligands |

### Switching predictors in `seq-fold-eval` and `protein-binder-design`

Set `predictor` to `boltz2` (default), `alphafold3`, `protenix`, or `chai1`. All four emit unified `confidence_score` and `iptm` columns so downstream `DESIGN_FILTER` expressions remain portable across oracles.

## 3. Additional Extended Modules

| Module | Best for | Key outputs & integration |
| --- | --- | --- |
| `local/modules/bindcraft` | One-shot AF2 hallucination + ProteinMPNN + filtering for protein-protein binders when you prefer a single-step target-hotspot trajectory over `RFD3 → MPNN` | `bindcraft_designs.tsv` (`complex_plddt`, `iptm`, `i_pae_normalized`, `dG`, `shape_complementarity`) + accepted PDBs |
| `local/modules/esmfold2` | High-throughput MSA-free single-sequence folding of monomers or rapid self-consistency pre-screening before Boltz-2 / AF3 | `esmfold2_metrics.tsv` (`plddt`, `complex_plddt`, `ptm`) directly consumable by `DESIGN_FILTER` |
| `local/modules/foldseek-easycluster` | Structural dereplication (`foldseek easy-cluster` → `cluster_id`, `cluster_representative`, `cluster_size`) | Auto-wired in `design_postprocess` (`run_foldseek_cluster: true`) or merge into `DESIGN_FILTER` with `group_by: cluster_id` and `max_per_group: 1` |
| `local/modules/foldseek-easysearch` | PDB/AFDB structural novelty screening (`foldseek easy-search` → `max_alntmscore`, `min_evalue`, `is_novel`) | Requires a reference `target_db`; merge `foldseek_search.tsv` into `DESIGN_FILTER` to gate on `is_novel == "true"` or `max_alntmscore < 0.5` |
