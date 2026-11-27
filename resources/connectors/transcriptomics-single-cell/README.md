Transcriptomics and Single-Cell unifies public bulk tissue expression, eQTL, and single-cell/spatial atlas metadata across the GTEx Portal, CZ CELLxGENE Discover, CZ CellGuide, the Human Cell Atlas (HCA), and HuBMAP.

### Capabilities

- Query GTEx reference genes, tissue sites, median tissue expression (in TPM), tissue-level eGenes, and single-tissue eQTL variants.
- Browse CZ CELLxGENE Discover collections and constituent datasets with cell counts, assays, and standardized UBERON, Cell Ontology (`CL`), and MONDO annotations.
- Look up canonical cell-type marker genes and supporting publications from the pinned CZ CellGuide snapshot by Cell Ontology ID.
- Search Human Cell Atlas projects via Azul by organ or project ID, inspecting species, specimen organs, estimated cell counts, and data-use restrictions.
- Search public HuBMAP datasets, samples, and donors across spatial and single-cell assay modalities.

### Research & Workflow Use Cases

- Establish baseline expression profiles across healthy human tissues in GTEx before comparing disease vs. control cohorts.
- Colocalize GWAS risk loci or candidate regulatory variants with GTEx single-tissue eQTLs and allelic fold-change estimates.
- Cross-check cluster annotations in scRNA-seq analyses against canonical CellGuide marker genes for matching `CL` ontology terms.
- Discover open reference single-cell and spatial datasets across CELLxGENE, HCA, and HuBMAP for label transfer, integration, or replication.

### Tools

| Tool | Description |
| --- | --- |
| `usage_terms` | Return usage policies and data access notes for GTEx, CZ CELLxGENE, CellGuide, Human Cell Atlas, and HuBMAP. |
| `gtex_search_genes` | Search GTEx reference genes by gene symbol or Ensembl ID. |
| `gtex_list_tissues` | List GTEx tissue sites with RNA-seq sample counts and UBERON ontology mappings. |
| `gtex_gene_expression` | Fetch median GTEx gene expression (in TPM) across tissues or for a specific tissue site. |
| `gtex_list_egenes` | List GTEx eGenes and log2 allelic fold changes for a specified tissue site. |
| `gtex_list_eqtls` | List single-tissue GTEx eQTL variants, normalized effect sizes, and p-values for a gene and tissue. |
| `cellxgene_list_collections` | List public CZ CELLxGENE Discover collections with DOIs and dataset counts. |
| `cellxgene_get_collection` | Fetch metadata for a single CELLxGENE collection by collection UUID. |
| `cellxgene_list_datasets` | List single-cell datasets embedded in a CELLxGENE collection, including cell counts, assays, and ontology terms. |
| `cellxgene_get_dataset` | Fetch metadata and Explorer URL for a single dataset within a CELLxGENE collection. |
| `cellguide_cell_type` | Look up canonical marker genes and supporting publications for a Cell Ontology (`CL`) ID in the pinned CellGuide snapshot. |
| `hca_search_projects` | Search Human Cell Atlas projects via Azul, optionally filtered by organ or project ID. |
| `hca_get_project` | Fetch metadata, species, organ annotations, and access status for a single HCA project. |
| `hubmap_search_datasets` | Search public HuBMAP datasets by keyword, returning assay types and portal links. |
| `hubmap_search_samples` | Search HuBMAP sample metadata and linked donor identifiers. |
| `hubmap_search_donors` | Search HuBMAP donor metadata and access levels. |

### Usage Notes

- No sign-in or API key is required. HCA searches via Azul support `offset + limit <= 100`, and CellGuide queries read from a pinned public snapshot (`CELLGUIDE_SNAPSHOT`).
- Read-only aggregate and metadata access only: the server does not download `h5ad` matrices or raw sequencing files, and controlled-access HCA (non-`NRES`) or non-public HuBMAP records omit data URLs.
- Homepage: [Phi Packages](https://github.com/dxsbiocc/phi-packages). Cite [GTEx](https://gtexportal.org/), [CZ CELLxGENE](https://cellxgene.cziscience.com/), [Human Cell Atlas](https://data.humancellatlas.org/), and [HuBMAP](https://portal.hubmapconsortium.org/) alongside the primary study authors.
