# Campaign Scaling & Iterative Refinement Guide

## 1. Parallel Chunking (`n_chunks` vs `n_batches`)

In `local/workflows/protein-binder-design`, total generated backbones equal:

`total_backbones = n_chunks × n_batches × diffusion_batch_size`

- **Smoke test**: `n_chunks: 1`, `n_batches: 1`, `diffusion_batch_size: 2`, `num_seq_per_pdb: 2` (4 total designs).
- **Pilot**: `n_chunks: 2`, `n_batches: 2`, `diffusion_batch_size: 8`, `num_seq_per_pdb: 4` (128 total designs).
- **Campaign**: Increase `n_chunks` to fan out independent GPU jobs across a cluster (e.g. `n_chunks: 10`, `n_batches: 5`, `diffusion_batch_size: 8` → 400 backbones × 8 sequences = 3,200 designs), while keeping each chunk small enough to fit walltime and GPU memory limits.
- All chunk metric tables (`sequences.tsv`, predictor `metrics.tsv`, and `protein_metrics.tsv`) automatically converge into a single campaign-level `DESIGN_FILTER` call at the end.

## 2. Agent-Driven Iterative Refinement (`engine.cycle()` Equivalent)

Because Nextflow DSL2 executes a static DAG, run iterative MPNN ↔ Oracle optimization as an agent-managed outer loop:

1. **Round 1**: Run `local/workflows/protein-binder-design` (or `local/subworkflows/seq-fold-eval` + `local/modules/design-filter`) with exploratory MPNN temperature (`0.2`).
2. **Diagnose**: Run `pd_funnel_report` on Round 1 outputs to identify the top-performing predicted structures and failure bottlenecks.
3. **Round 2+**: Pass the Round 1 passing predicted structures (`.cif`/`.pdb`) as `structures` into `local/subworkflows/seq-fold-eval` with lower temperature (`0.1`) and optional `fixed_residues` locking high-confidence interface contacts, then re-evaluate with `local/modules/design-filter`.
4. **Stop rule**: Stop when top-N `iptm` / `complex_plddt` and `ca_rmsd` stop improving across rounds or meet the campaign acceptance criteria.
