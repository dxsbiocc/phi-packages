---
name: structure-input-prep
description: >-
  Prepare and validate protein/complex structure inputs before launching GPU
  design wrappers: inspect PDB/mmCIF chain IDs, residue numbering gaps, and
  HETATM ligands with `pd_inspect_structure`; validate RFdiffusion3 JSON/contig
  specifications against target structures with `pd_validate_rfd3_spec`; and
  follow the contig, chain-break, and ligand CCD/SMILES authoring guide. Use
  whenever preparing a target PDB/CIF, authoring RFD3 contigs, or debugging
  chain/ligand parsing errors.
phi:
  environment: phi:python@1
  attachTo:
    - ProteinDesign
  scripts:
    - name: inspect_structure
      description: >-
        Inspect a PDB or mmCIF structure file and report exact protein chain
        IDs, residue numbering ranges, internal numbering gaps, polymer length,
        and non-water HETATM ligands (name, chain, atom count). Run before
        authoring any RFD3 contig or chain selection.
      run: [python3, ./scripts/inspect_structure.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: Path to a PDB (.pdb) or mmCIF (.cif/.mmcif) structure file.
      approval: read
      output: ./schemas/inspect-structure-result.json
      timeoutSeconds: 60
    - name: validate_rfd3_spec
      description: >-
        Validate an RFdiffusion3 input JSON file or a quick contig string
        (optionally cross-checking referenced chains, residue ranges, and
        ligand codes against an input PDB/mmCIF structure) before submitting a
        GPU wrapper run.
      run: [python3, ./scripts/validate_rfd3_spec.py]
      args:
        type: object
        additionalProperties: false
        properties:
          inputs_json:
            type: string
            format: input-path
            description: Optional path to an RFD3 JSON specification file.
          contig:
            type: string
            description: Optional RFD3 contig string to validate directly (e.g. "A1-100/0 60-90").
          input_structure:
            type: string
            format: input-path
            description: Optional target PDB/mmCIF structure to cross-check chain IDs, residue ranges, and ligands.
          ligand:
            type: string
            description: Optional comma-separated ligand residue code(s) to verify against input_structure.
      approval: read
      output: ./schemas/validate-rfd3-spec-result.json
      timeoutSeconds: 60
---

# Structure Input Preparation & Validation

Always run deterministic input validation before launching `local/modules/rfd3`,
`local/modules/ligandmpnn`, `local/modules/boltz2`, or `local/modules/alphafold3`.

## Reference guide

- [references/contig-and-ligand-guide.md](references/contig-and-ligand-guide.md) — RFD3 contig syntax, chain breaks (`/0`), motif scaffolding, ligand CCD vs SMILES rules, and common input troubleshooting.

## Standard prep workflow

1. Call `pd_inspect_structure` on the target PDB/mmCIF to confirm real chain IDs, residue numbering ranges, missing-density gaps, and HETATM ligand codes.
2. Author the RFD3 `contig` or `inputs_json` using the exact chain IDs and residue numbers reported by `pd_inspect_structure`.
3. Call `pd_validate_rfd3_spec` to verify syntax and cross-check every referenced chain, residue range, and ligand code against the structure file.
4. Only after `valid: true` (zero errors), pass the validated inputs to `wrapper_plan_run`.
