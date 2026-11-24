# Contig & Ligand Preparation Guide

## 1. RFdiffusion3 Contig Syntax

Contigs describe which residues come from an input structure and which segments are generated de novo:

- **Fixed target chain segment**: `A1-120` (takes chain A residues 1 through 120 from `input_structure`).
- **De novo designed segment**: `60-90` (samples a new chain or linker between 60 and 90 residues long) or `75` (exact length 75).
- **Chain break**: `/0` separates independent polypeptide chains with no peptide bond between them.
- **Binder against fixed target**: `A1-120/0 60-85` keeps target chain A (residues 1–120) fixed and diffuses a separate 60–85 aa binder chain.
- **Motif scaffolding**: `20-35/A40-60/15-25` embeds fixed motif `A40-60` inside a newly diffused scaffold.

## 2. Common Input Pitfalls

1. **Missing residues / numbering gaps**: Crystal structures often have unmodeled loops (e.g. `A1-45` then `A52-120`). Run `pd_inspect_structure` first so you know where gaps exist instead of assuming `A1-120` is contiguous.
2. **Non-standard chain IDs**: Downloaded PDB files may use chain `B` or `C` instead of `A`. Never hardcode `A` without inspecting.
3. **Multiple conformers / altlocs**: Ensure you know whether alternative conformations exist; our metric scripts keep default conformer (`""`/`"A"`).
4. **Water and buffer ions**: `HOH` waters are ignored by `pd_inspect_structure`, while buffer ligands (`SO4`, `GOL`, `EDO`) appear under `ligands` — verify you pass only the intended biological ligand code to `ligand`.

## 3. Ligand Specification Across Tools

- **RFD3**: Pass the exact 3-letter HETATM residue name from `input_structure` via `ligand` (e.g. `"ATP"` or `"HEM"`).
- **LigandMPNN**: Reads the HETATM ligand directly from the input PDB/CIF when `model_type: ligand_mpnn`.
- **Boltz-2 & AlphaFold 3**: Prefer official **CCD codes** (`ligand_ccd: "ATP"`) over SMILES whenever the ligand exists in the PDB Chemical Component Dictionary; fall back to `ligand_smiles` for novel small molecules not in CCD.
