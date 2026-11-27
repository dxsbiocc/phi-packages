STRING is a database of known and predicted protein-protein interactions covering thousands of organisms, integrating physical binding experiments, curated pathway databases, co-expression, genomic neighborhood, and literature text mining.

### Capabilities

- Resolve gene symbols or protein identifiers within a target species and retrieve interaction partners and local network neighborhoods.
- Inspect combined interaction confidence scores alongside channel-specific scores (experiments, curated databases, co-expression, text mining, co-occurrence, neighborhood, fusion).
- Query cross-species homologs and orthologous groups to map findings between model organisms and humans.
- Run functional enrichment analysis on protein sets across GO, KEGG, Reactome, WikiPathways, UniProt keywords, and disease ontologies.
- Export node and edge tables suitable for downstream network visualization in tools such as Cytoscape or igraph.

### Research & Workflow Use Cases

- Build physical and functional interaction subnetworks from differentially expressed genes or proteomics hits to identify hub proteins and complexes.
- Distinguish experimentally validated physical interactions from text-mining or co-expression associations before designing co-IP or proximity-labeling assays.
- Check whether candidate interaction modules are conserved across mouse, zebrafish, fly, and human orthologs.
- Interpret gene modules or clustering results via pathway and domain enrichment statistics.

### Usage Notes

- No sign-in or API key is required to query the official STRING MCP endpoint.
- Interaction scores and enrichment backgrounds depend on the active STRING version and confidence threshold; record both parameters for reproducibility.
- Homepage: [STRING MCP](https://github.com/meringlab/string-mcp). Predicted functional associations do not imply direct physical binding; verify critical edges against experimental evidence and cite [STRING](https://string-db.org/).
