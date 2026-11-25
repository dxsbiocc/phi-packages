# Changelog — protein-design

## 0.1.0 (2026-11-23)

Initial release: end-to-end de novo protein and small-molecule binder design plugin (`toolPrefix: pd`, agent `ProteinDesign`, shared `phi:python@1` environment).

- `protein-design-workflow`: playbook, block selection, scaling, and campaign iteration guides paired with local Nextflow wrappers (`rfd3`, `ligandmpnn`, `boltz2`, `alphafold3`, `protenix`, `chai1`, `esmfold2`, `bindcraft`, `foldseek`, `protein_metrics`, `design_filter`, `seq_fold_eval`, `design_postprocess`, `protein_binder_design`).
- `structure-input-prep`: target inspection (`pd_inspect_structure`) and RFdiffusion3 contig/ligand specification validator (`pd_validate_rfd3_spec`).
- `design-evaluation`: confidence, self-consistency RMSD, interface, and affinity funnel diagnostics (`pd_funnel_report`).
