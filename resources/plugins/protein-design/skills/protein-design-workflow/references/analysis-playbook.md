# Protein Design Analysis Playbook

Use this playbook to advance any protein design campaign through four disciplined stages.

## Stage 1 — Target & Constraint Inspection

- **Call**: `pd_inspect_structure` on the target PDB/mmCIF, then `pd_validate_rfd3_spec` on your candidate contig or JSON specification.
- **Read**: Exact chain IDs, continuous residue ranges (watch for missing density gaps), and ligand HETATM residue names.
- **Decide**: Choose `contig` length (typically 55–100 aa for compact de novo PPI binders; 90–150 aa for small-molecule pocket binders) and whether to fix catalytic/interface motifs.
- **Red lines**:
  - Never guess chain IDs or residue numbers from memory or a paper figure — always inspect the actual PDB/mmCIF file.
  - Never pass water (`HOH`) or crystallization buffer ions as target ligands unless intentionally designing a metal-binding site.

## Stage 2 — Smoke Test (End-to-End Plumbing)

- **Call**: `local/workflows/protein-binder-design` with `n_chunks: 1`, `n_batches: 1`, `diffusion_batch_size: 2`, `num_seq_per_pdb: 2`, `recycling_steps: 3`, `diffusion_samples: 1`, `filter_expr: ""`.
- **Read**: Confirm `01_rfd3`, `02_ligandmpnn`, `03_boltz2` (or `03_alphafold3`), `04_metrics`, and `05_selected` all produce non-empty outputs and that `design_id` / `parent_design_id` join cleanly in `merged.tsv`.
- **Decide**: Fix any chain-order mismatch, missing `target_sequence`, or ligand CCD/SMILES mismatch before running a pilot.
- **Red lines**:
  - A smoke test with 0 passing designs under strict filters is expected — leave `filter_expr: ""` during the smoke test to inspect raw metric distributions.

## Stage 3 — Pilot & Calibration

- **Call**: Run 16–32 backbones × 4–8 sequences (`n_chunks: 2`, `n_batches: 2`, `diffusion_batch_size: 8`, `num_seq_per_pdb: 4`) and inspect results with `pd_funnel_report`.
- **Read**: Pass rate, median `complex_plddt` / `plddt`, `iptm`, `min_pae_interaction`, `ca_rmsd`, `interface_contacts`, and `clash_count`.
- **Decide**:
  - High `ca_rmsd` (> 2.5 Å) but good `plddt`: adjust MPNN `temperature` (lower to `0.1`) or increase `num_seq_per_pdb`.
  - Low `iptm` (< 0.6) or few `interface_contacts`: adjust binder length/contig or add hotspot constraints in `inputs_json`.
  - High `clash_count`: ensure `use_potentials: true` in Boltz-2 or inspect ligand geometry.
- **Red lines**:
  - Do not scale a zero-pass pilot branch to production without diagnosing and fixing the root cause first.

## Stage 4 — Campaign Scale & Diverse Selection

- **Call**: Scale `n_chunks` and `n_batches` based on pilot pass rate to reach the target number of orderable designs (typically 24–96 diverse sequences).
- **Read**: Run `pd_funnel_report` on final `merged.tsv`, `filtered.tsv`, and `summary.json`.
- **Decide**: Enforce `group_by: parent_design_id` and `max_per_group: 1` or `2` in `DESIGN_FILTER` so final selections represent distinct backbones rather than near-duplicate sequences from one scaffold.
- **Red lines**:
  - Final ranking must compare metrics measured by the **same** structure predictor (`boltz2` vs `alphafold3` scores are not directly interchangeable).
