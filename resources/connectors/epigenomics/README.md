Epigenomics queries the ENCODE Portal REST API to surface public experiment, file, and biosample metadata from the ENCODE Consortium and the NIH Roadmap Epigenomics Program.

### Capabilities

- Search released ENCODE experiments by free-text keyword or assay title (such as `ChIP-seq`, `ATAC-seq`, `DNase-seq`, or `WGBS`).
- Filter specifically for Roadmap Epigenomics reference datasets (`award.rfa=Roadmap`) across human tissues and primary cell types.
- Query file metadata by keyword, file format, or `ENCFF` accession, returning direct download URLs for released unrestricted files.
- Search and inspect `ENCBS` biosample records with UBERON tissue and Cell Ontology (`CL`) annotations, donor links, and species identifiers.

### Research & Workflow Use Cases

- Locate public transcription factor ChIP-seq, histone modification, and chromatin accessibility tracks for regulatory element annotation.
- Match reference epigenomes by tissue, primary cell ontology term, and assay before running peak colocalization or motif enrichment pipelines.
- Audit file formats (`bigWig`, `bed`, `bam`), output types, and release status prior to batch downloading from ENCODE.
- Compare chromatin state and accessibility datasets across human organs profiled in Roadmap Epigenomics.

### Tools

| Tool | Description |
| --- | --- |
| `usage_terms` | Return data usage policies and source notes for ENCODE and Roadmap Epigenomics metadata. |
| `encode_search_experiments` | Search released ENCODE experiments by keyword or assay title (such as `ChIP-seq` or `ATAC-seq`). |
| `encode_get_experiment` | Fetch metadata for a single ENCODE experiment by accession (for example, `ENCSR818VNL`). |
| `encode_search_roadmap` | Search released ENCODE experiments funded under the Roadmap Epigenomics program (`award.rfa=Roadmap`). |
| `encode_search_files` | Search released ENCODE file metadata by keyword or file format, omitting URLs for restricted files. |
| `encode_get_file` | Fetch metadata for a single ENCODE file by `ENCFF` accession. |
| `encode_search_biosamples` | Search released ENCODE biosamples with UBERON tissue and Cell Ontology annotations. |
| `encode_get_biosample` | Fetch metadata for a single released ENCODE biosample by `ENCBS` accession. |

### Usage Notes

- No sign-in or API key is required. Search pagination via the ENCODE portal supports `offset + limit <= 100`.
- Read-only metadata only: the server does not download raw file bytes, and restricted or unreleased ENCODE files return `url: null`.
- Homepage: [Phi Packages](https://github.com/dxsbiocc/phi-packages). Follow [ENCODE data use and citation guidelines](https://www.encodeproject.org/help/getting-started/) when publishing analyses.
