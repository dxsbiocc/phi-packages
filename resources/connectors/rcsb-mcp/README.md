RCSB PDB is the US data center for the global Protein Data Bank archive, providing open access to experimentally determined 3D structures of proteins, nucleic acids, macromolecular assemblies, and bound small-molecule ligands.

### Capabilities

- Search PDB entries by macromolecule name, UniProt/sequence query, chemical ligand ID, source organism, or experimental method (X-ray crystallography, cryo-EM, NMR).
- Inspect entry metadata including resolution, experimental conditions, polymer chains, deposition dates, and primary citations.
- Retrieve polymer sequences, ligand chemical descriptors, chain stoichiometries, and binding-site annotations.
- Compare multiple structures of the same protein or complex across conformations, mutants, and co-crystallized ligands.

### Research & Workflow Use Cases

- Map missense variants, post-translational modifications, and functional domains onto experimental 3D structures.
- Select high-resolution receptor templates and verify co-crystallized ligands before molecular docking, free-energy calculations, or molecular dynamics simulations.
- Compare active-site geometry and binding-pocket conservation across orthologs or protein family members.
- Compile curated PDB accession tables with resolution, ligand IDs, and literature references for structural biology or medicinal chemistry studies.

### Usage Notes

- No sign-in or API key is required to query the public RCSB PDB MCP server.
- This connector currently points to RCSB's beta MCP endpoint (`mcp-beta.rcsb.org`); tool schemas and response fields may evolve across upstream updates.
- Homepage: [RCSB PDB MCP](https://github.com/rcsb/rcsb-mcp). Record PDB entry IDs, structure revision versions, and primary citations per [RCSB PDB policies](https://www.rcsb.org/pages/usage-policy).
