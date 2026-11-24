# Block & Tool Selection Guide

## 1. Sequence Design (`local/modules/ligandmpnn`)

| `model_type` | Best for | Recommended `temperature` |
| --- | --- | --- |
| `protein_mpnn` | Standard protein-only backbones and PPI binders | `0.1` (production) / `0.2` (exploration) |
| `ligand_mpnn` | Small-molecule binders, enzymes, metal/cofactor complexes, nucleic-acid interfaces | `0.1`–`0.2` |
| `soluble_mpnn` | Improving solubility and expression of hydrophobic de novo scaffolds | `0.1`–`0.2` |

- Keep `output_mode: sequences` by default: downstream `boltz2` and `alphafold3` modules read `sequences.tsv` directly without needing per-sequence packed PDBs.
- Use `output_mode: structures` only when you explicitly need side-chain-packed PDBs for physics scoring or visual inspection before refolding.

## 2. Structure Prediction (`boltz2` vs `alphafold3`)

| Dimension | `local/modules/boltz2` | `local/modules/alphafold3` |
| --- | --- | --- |
| License & weights | MIT code & open weights (`boltz_cache_dir`) | Apache-2.0 code; **non-commercial** weights supplied by user (`af3_model_dir`) |
| Ligand affinity prediction | Supported (`predict_affinity: true` -> `affinity_pred_value`, `affinity_probability_binary`) | Not supported (structure & confidence only) |
| Fast screening mode | `use_msa_server: false`, `recycling_steps: 3`, `diffusion_samples: 1`, `use_potentials: true` | `run_data_pipeline: false`, `model_seeds: "1"`, `num_diffusion_samples: 5` |
| Key confidence columns | `confidence_score`, `complex_plddt`, `iptm`, `ligand_iptm`, `complex_ipde` | `confidence_score` (`ranking_score`), `complex_plddt` (`plddt/100`), `iptm`, `min_pae_interaction`, `has_clash` |

### Switching predictors in `seq-fold-eval` and `protein-binder-design`

Set `predictor: boltz2` (default) for fast open-license screening and small-molecule affinity estimation, or set `predictor: alphafold3` with `af3_model_dir` for final orthogonal validation on top candidates. Both predictors emit `confidence_score`, `complex_plddt`, and `iptm` so downstream `DESIGN_FILTER` expressions remain portable.

## 3. Additional Extended Modules

| Module | Best for | Key outputs & integration |
| --- | --- | --- |
| `local/modules/bindcraft` | One-shot AF2 hallucination + ProteinMPNN + filtering for protein-protein binders when you prefer a single-step target-hotspot trajectory over `RFD3 → MPNN` | `bindcraft_designs.tsv` (`complex_plddt`, `iptm`, `i_pae_normalized`, `dG`, `shape_complementarity`) + accepted PDBs |
| `local/modules/esmfold2` | High-throughput MSA-free single-sequence folding of monomers or rapid self-consistency pre-screening before Boltz-2 / AF3 | `esmfold2_metrics.tsv` (`plddt`, `complex_plddt`, `ptm`) directly consumable by `DESIGN_FILTER` |
| `local/modules/foldseek` | Structural dereplication (`mode: cluster` → `cluster_id`) and PDB/AFDB structural novelty screening (`mode: search` → `max_alntmscore`, `is_novel`) | Merge `foldseek.tsv` into `DESIGN_FILTER` with `group_by: cluster_id` and `max_per_group: 1` for structural diversity |
