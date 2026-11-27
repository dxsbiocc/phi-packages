# Changelog — protein-design

## 0.2.0 (2026-11-26)

- `design-evaluation`: add pure-stdlib `pd_developability_audit` (`scripts/developability_audit.py` + `schemas/developability-audit-result.json`) for pre-wetlab biophysical & chemical-liability triage on FASTA or `DESIGN_FILTER` tables (average/monoisotopic MW, Henderson-Hasselbalch $pI$, net charge & $|pI - \text{pH}| < 1.0$ precipitation warning, reduced/oxidized $\varepsilon_{280}$ & $A_{280}^{0.1\%}$ with UV280-blind alert, GRAVY, Guruprasad Instability Index, Aliphatic Index, 9-mer max Kyte-Doolittle hydrophobic patch, unpaired Cys, `NG/NS/NT` deamidation, `DG/DS` isomerization, `N-X-[ST]` N-glycosylation sequons, N-terminal processing notes, and `pass/caution/high_risk` tiers).
- `agents/ProteinDesign.md` & `references/metrics-and-thresholds.md`: add Pre-Wetlab Critic & Developability Gate before ordering DNA.

## 0.1.0 (2026-11-23)

Initial release: end-to-end de novo protein and small-molecule binder design plugin (`toolPrefix: pd`, agent `ProteinDesign`, shared `phi:python@1` environment).

- `protein-design-workflow`: playbook, block selection, scaling, and campaign iteration guides paired with local Nextflow wrappers (`rfd3`, `ligandmpnn`, `boltz2`, `alphafold3`, `protenix`, `chai1`, `esmfold2`, `bindcraft`, `foldseek`, `protein_metrics`, `design_filter`, `seq_fold_eval`, `design_postprocess`, `protein_binder_design`).
- `structure-input-prep`: target inspection (`pd_inspect_structure`) and RFdiffusion3 contig/ligand specification validator (`pd_validate_rfd3_spec`).
- `design-evaluation`: confidence, self-consistency RMSD, interface, and affinity funnel diagnostics (`pd_funnel_report`).
