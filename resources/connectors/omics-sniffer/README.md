Omics Sniffer is a cross-repository discovery engine indexing standardized dataset, sample, and file artifact metadata across 30+ public life-science archives, including GEO, SRA, ArrayExpress, BioStudies, PRIDE, MassIVE, CZ CELLxGENE, ENCODE, GTEx, GDC/TCGA, cBioPortal, CPTAC, DepMap, HuBMAP, Human Protein Atlas, ChIP-Atlas, GWAS Catalog, FinnGen, MetaboLights, Metabolomics Workbench, Zenodo, and Figshare.

### Capabilities

- Inspect active Parquet catalog release metadata, total indexed dataset counts, per-repository coverage, and incremental sync timestamps.
- Filter public datasets simultaneously by keyword, organism, assay/modality, disease, tissue/organ, source repository, access tier, and minimum sample count.
- Look up complete dataset records by canonical ID (e.g., `geo:GSE164760`) or bare repository accession (`GSE...`, `E-MTAB-...`, `PXD...`, `TCGA-...`), including cross-repository member accessions and linked PubMed/DOI publications.
- Browse paginated sample-level attributes (organism, tissue, disease, experimental condition, sex, age, and SRA run IDs) and downloadable raw or processed file artifacts (`counts`, `TPM`, `h5ad`, `mtx`, `FASTQ`/FTP roots, clinical tables).

### Research & Workflow Use Cases

- Survey all publicly available bulk RNA-seq, single-cell, spatial, epigenomic, proteomic, and metabolomic cohorts for a disease, organ, or rare condition in one query.
- Pre-screen cohorts that meet statistical power thresholds (such as `min_samples >= 50`) and provide open processed count or `h5ad` matrices ready for downstream bioinformatics workflows.
- Resolve cross-linked accessions (GEO, SRA, BioProject, PRIDE) and deposited supplementary matrices associated with a published study.

### Tools

| Tool | Description |
| --- | --- |
| `usage_terms` | Return usage terms, repository coverage, and active Parquet snapshot status for Omics Sniffer. |
| `get_catalog_stats` | Return active Parquet snapshot release info, total indexed datasets, per-repository dataset counts, and incremental auto-sync freshness. |
| `search_omics_datasets` | Search curated public multi-omics datasets across 30+ repositories by keyword, organism, assay, disease, tissue, source, minimum sample count, and access level. |
| `get_dataset` | Fetch full metadata, cross-repository member accessions, and linked PubMed/DOI publications for one dataset by canonical ID or bare accession. |
| `list_samples` | List paginated sample-level metadata (organism, tissue, disease, condition, sex, age, attributes, and SRA run IDs) for a dataset. |
| `list_artifacts` | List raw, processed, and supplementary file artifacts (`counts`, `TPM`, `h5ad`, `FASTQ`/FTP roots, clinical tables) for a dataset, omitting download URLs for controlled-access artifacts. |

### Usage Notes

- No sign-in or API key is required.
- Queries run against read-only Parquet catalog snapshots updated incrementally; controlled-access records return metadata only and omit download URLs.
- Homepage: [Phi Packages](https://github.com/dxsbiocc/phi-packages). Comply with the data use licenses and citation policies of each originating repository and study author when reusing datasets.
