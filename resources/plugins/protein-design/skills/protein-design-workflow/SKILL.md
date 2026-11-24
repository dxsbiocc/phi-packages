---
name: protein-design-workflow
description: >-
  End-to-end tutorial and decision playbook for computational protein design in
  Phi: choosing between de novo scaffold generation, PPI binder design,
  small-molecule binder design, and sequence redesign; mapping tasks to the
  local Nextflow wrappers (rfd3, ligandmpnn, boltz2, alphafold3,
  protein-metrics, design-filter, seq-fold-eval, design-postprocess,
  protein-binder-design); scaling from smoke tests to production campaigns; and
  running agent-driven iterative refinement cycles. Use whenever planning a
  protein design pipeline, selecting tools/predictors, or explaining how
  Chaperone-style block workflows map onto Phi Nextflow wrappers.
phi:
  attachTo:
    - ProteinDesign
---

# Protein Design Workflow Guide

This skill provides the end-to-end decision playbook for running protein and
biomolecular binder design campaigns in Phi using the local Nextflow wrappers.

## Reference guides

1. [references/analysis-playbook.md](references/analysis-playbook.md) — Stage-by-stage Call → Read → Decide → Red Lines procedure from target inspection to diverse candidate handoff.
2. [references/block-selection-guide.md](references/block-selection-guide.md) — How to choose backbone generators, MPNN variants (`protein_mpnn`, `ligand_mpnn`, `soluble_mpnn`), and structure prediction oracles (`boltz2` vs `alphafold3`).
3. [references/campaign-scaling-guide.md](references/campaign-scaling-guide.md) — Moving from a 1–2 design smoke test to chunked production campaigns (`n_chunks`, `n_batches`, `output_mode: sequences`) and agent-driven iterative refinement cycles.

## Canonical wrapper map

| Design goal | Recommended Wrapper | Key inputs |
| --- | --- | --- |
| End-to-end binder campaign (`RFD3 → MPNN → Boltz2/AF3 → Metrics → Filter`) | `local/workflows/protein-binder-design` | `input_structure` + `contig` (or `inputs_json`), `target_sequence`, weight dirs |
| Sequence redesign + refolding + scoring on existing backbones | `local/subworkflows/seq-fold-eval` | `structures`, `predictor` (`boltz2` or `alphafold3`), `num_seq_per_pdb` |
| Score & filter already-folded structures + upstream metric tables | `local/subworkflows/design-postprocess` | `structures`, `extra_tables`, `expr`, `rank_by`, `group_by` |
| Single-stage backbone generation | `local/modules/rfd3` | `inputs_json` or `contig` + `input_structure` |
| Single-stage inverse folding | `local/modules/ligandmpnn` | `structures`, `model_type`, `temperature` |
| Single-stage structure prediction | `local/modules/boltz2` or `local/modules/alphafold3` | `inputs` or `sequences_tsv`, cache/model dir |

## Hard rules

- Always validate inputs with `pd_inspect_structure` and `pd_validate_rfd3_spec` (from `structure-input-prep`) before submitting GPU wrappers.
- Never scale straight to 100+ designs without a passing smoke test (`n_chunks=1`, `n_batches=1`, `diffusion_batch_size=2`, `num_seq_per_pdb=2`).
- Keep `mpnn_output_mode: sequences` for campaigns so intermediate per-sequence PDB files do not flood disk.
- Confirm model weights and license constraints (`af3_model_dir` is non-commercial and user-supplied) before launching any GPU wrapper.
