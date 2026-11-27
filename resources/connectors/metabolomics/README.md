Metabolomics connects to the EMBL-EBI MetaboLights repository to retrieve public metabolomics study accessions, ISA-Tab study descriptions, release status flags, assay metadata files, and public FTP file manifests.

### Capabilities

- List public MetaboLights study accessions (`MTBLS...`) and inspect individual study titles, descriptions, and release status.
- List files associated with a study, including file names, types, and public EMBL-EBI FTP download URLs.
- List ISA-Tab assay files (`a_*.txt`) within a study to inspect analytical modalities and experimental runs.
- Distinguish public studies from private submissions, omitting download URLs when a study is marked private.

### Research & Workflow Use Cases

- Locate open LC-MS, GC-MS, and NMR metabolomics datasets and study designs for reanalysis or benchmarking.
- Review assay layouts and file inventories from comparable studies before designing sample-prep or mass spectrometry protocols.
- Assemble public reference cohorts for biomarker discovery, exposomics, or metabolic pathway studies.
- Verify MetaboLights accession numbers and deposited file manifests cited in manuscripts.

### Tools

| Tool | Description |
| --- | --- |
| `usage_terms` | Return data usage policies and API notes for EMBL-EBI MetaboLights. |
| `metabolights_list_studies` | List public MetaboLights study accessions (`MTBLS...`) with pagination. |
| `metabolights_get_study` | Fetch the title, description, and release status for a single MetaboLights study by `MTBLS` accession. |
| `metabolights_list_files` | List file records and public FTP URLs for a MetaboLights study, omitting URLs for private studies. |
| `metabolights_list_assays` | List assay metadata file records associated with a MetaboLights study. |

### Usage Notes

- No sign-in or API key is required.
- Read-only metadata only: the server does not download raw spectral files or expose parsed metabolite abundance matrices (metabolite assignment tables `m_*.tsv` must be retrieved from the listed study files).
- Homepage: [Phi Packages](https://github.com/dxsbiocc/phi-packages). Follow [EMBL-EBI MetaboLights](https://www.ebi.ac.uk/metabolights/about) data reuse and citation guidelines.
