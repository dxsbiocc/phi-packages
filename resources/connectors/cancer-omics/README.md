Cancer Omics aggregates open cancer genomics and proteomics metadata from the NCI Genomic Data Commons (GDC), the public cBioPortal REST API, and the NCI Proteomic Data Commons (PDC), returning normalized donor, sample, gene, variant, dataset, and file records through a single MCP endpoint.

### Capabilities

- Search and inspect NCI GDC donor cases by case UUID or submitter ID, including disease type, primary site, and age at index.
- Query GDC gene annotations (GRCh38), simple somatic mutations (SSMs), and file metadata with explicit open vs. controlled access flags.
- Browse public cBioPortal studies, patient cohorts, sample types, and HUGO/Entrez gene records.
- List and inspect PDC proteogenomic studies across CPTAC and related programs, including disease types and analytical fractions.

### Research & Workflow Use Cases

- Screen TCGA and GDC cohorts to verify sample counts, data formats, and open vs. controlled file availability before downloading data.
- Compare study designs, sample sizes, and reference genomes across public cBioPortal cohorts during target prioritization.
- Pair CPTAC proteomics studies in PDC with matching genomic or transcriptomic cohorts for multi-omics integration.
- Assemble harmonized cross-repository study and case manifests for survival, differential expression, or somatic mutation analyses.

### Tools

| Tool | Description |
| --- | --- |
| `usage_terms` | Return data usage policies and source links for NCI GDC, cBioPortal, and PDC. |
| `gdc_search_cases` | Search NCI GDC donor cases by case UUID or submitter ID, returning disease type, primary site, and age at index. |
| `gdc_get_case` | Fetch a single GDC donor case record by case UUID or submitter ID. |
| `gdc_search_files` | Search GDC file metadata filtered by case ID or access tier (`open` or `controlled`). |
| `gdc_get_gene` | Fetch a GDC gene record and GRCh38 genomic coordinates by gene symbol or Ensembl ID. |
| `gdc_search_mutations` | Search GDC simple somatic mutation (SSM) variant records, optionally filtered by Ensembl gene ID. |
| `cbioportal_list_studies` | List public cBioPortal cancer genomics studies with sample counts, cancer types, and reference genomes. |
| `cbioportal_get_study` | Fetch metadata for a single cBioPortal study by `study_id`. |
| `cbioportal_list_samples` | List sample identifiers and sample types within a cBioPortal study. |
| `cbioportal_list_patients` | List patient donor records associated with a cBioPortal study. |
| `cbioportal_get_gene` | Fetch a cBioPortal gene record and NCBI Entrez Gene CURIE by HUGO symbol. |
| `pdc_list_studies` | List NCI Proteomic Data Commons studies across CPTAC and related proteogenomic programs. |
| `pdc_get_study` | Fetch metadata for a single PDC study by `PDC` accession (for example, `PDC000759`). |

### Usage Notes

- No sign-in or API key is required. This connector queries the public cBioPortal REST API and is independent of the standalone cBioPortal OAuth connector.
- Read-only metadata only: the server does not download raw file bytes, and controlled-access GDC or cBioPortal records return `url: null`.
- Homepage: [Phi Packages](https://github.com/dxsbiocc/phi-packages). Follow the data access policies and citation guidelines of [NCI GDC](https://gdc.cancer.gov/), [cBioPortal](https://www.cbioportal.org/), and [NCI PDC](https://pdc.cancer.gov/).
