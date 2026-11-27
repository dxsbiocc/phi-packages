Proteomics connects to the EMBL-EBI PRIDE Archive and the Human Protein Atlas (HPA), combining public mass-spectrometry dataset and file discovery with human gene tissue-specificity and UniProt protein annotations.

### Capabilities

- Search public PRIDE Archive projects by keyword and inspect project titles, descriptions, submission dates, and `PXD` accessions.
- List file metadata within a PRIDE project, returning HTTPS download URLs converted from public PRIDE FTP locations.
- Search Human Protein Atlas genes by symbol or keyword, returning Ensembl gene CURIEs, GRCh38 coordinates, UniProt cross-references, and RNA tissue specificity labels.
- Retrieve single-gene HPA summaries and linked UniProt protein entries with protein class annotations by Ensembl ID.

### Research & Workflow Use Cases

- Discover reusable public mass-spectrometry datasets in PRIDE across specific tissues, diseases, or instrumentation workflows.
- Check baseline human RNA tissue specificity and distribution in HPA when evaluating candidate biomarkers or therapeutic targets.
- Cross-reference PRIDE project cohorts with HPA protein classes and canonical Ensembl/UniProt identifiers before downstream network or pathway analysis.

### Tools

| Tool | Description |
| --- | --- |
| `usage_terms` | Return data usage and licensing policies for PRIDE Archive and the Human Protein Atlas. |
| `pride_search_projects` | Search public PRIDE Archive proteomics projects by keyword. |
| `pride_get_project` | Fetch metadata for a single PRIDE project by `PXD` accession (for example, `PXD001357`). |
| `pride_list_files` | List file metadata and public HTTPS download URLs for a PRIDE project. |
| `hpa_search_genes` | Search Human Protein Atlas genes by keyword or symbol, returning Ensembl IDs, UniProt cross-references, and RNA tissue specificity. |
| `hpa_get_gene` | Fetch a single Human Protein Atlas gene summary and GRCh38 locus by Ensembl gene ID. |
| `hpa_list_proteins` | List UniProt protein records and protein classes linked to an Ensembl gene in the Human Protein Atlas. |

### Usage Notes

- No sign-in or API key is required.
- Read-only metadata only: the server does not download raw mass spectra, and non-public PRIDE files return `url: null`.
- Homepage: [Phi Packages](https://github.com/dxsbiocc/phi-packages). Cite the originating study alongside [PRIDE Archive](https://www.ebi.ac.uk/pride/) and the [Human Protein Atlas](https://www.proteinatlas.org/about/licence).
