OLS4 (Ontology Lookup Service v4) from EMBL-EBI indexes hundreds of biomedical and life-science ontologies, including Gene Ontology (GO), Cell Ontology (CL), Uberon, Mondo Disease Ontology (MONDO), Human Phenotype Ontology (HPO), EFO, and ChEBI.

### Capabilities

- Search ontology terms across all indexed ontologies or within a specific ontology, returning labels, definitions, synonyms, and CURIE/IRI identifiers.
- Traverse class hierarchies to inspect parent, child, ancestor, and descendant relationships.
- Find semantically similar classes across ontologies to assist with cross-ontology term alignment and harmonization.
- Map free-text cell type, tissue, phenotype, or disease descriptions to standardized ontology identifiers.
- Inspect ontology-level metadata, version strings, and release timestamps.

### Research & Workflow Use Cases

- Standardize cell type (`CL`), anatomical tissue (`UBERON`), and disease (`MONDO` / `EFO`) labels for single-cell atlases and multi-cohort metadata tables.
- Verify term definitions, exact synonyms, and parent-child relationships during pathway enrichment or knowledge-graph construction.
- Harmonize author-supplied sample descriptions from public studies into controlled vocabularies prior to cross-dataset integration.
- Check term obsolescence and replacement pointers when upgrading ontology versions in production pipelines.

### Usage Notes

- No sign-in or API key is required to query the public EMBL-EBI OLS4 MCP endpoint.
- Individual ontologies follow distinct release schedules and licenses; record the ontology version used in formal analyses.
- Homepage: [EMBL-EBI OLS4 MCP](https://www.ebi.ac.uk/ols4/mcp). Refer to OLS4 and each upstream ontology for licensing and citation details.
