cBioPortal for Cancer Genomics integrates multidimensional cancer genomics datasets from TCGA, ICGC, AACR Project GENIE, and institutional cohorts, exposing study metadata, patient and sample clinical attributes, and gene-level molecular alterations for cohort exploration.

### Capabilities

- Search published cancer genomics studies and inspect cancer types, sample counts, reference genomes, and linked publications.
- Browse patients and samples within a study to review sample types, clinical attributes, and cohort composition.
- Query somatic mutations and molecular profiles by gene across one or more tumor cohorts.
- Retrieve gene identifiers and annotations to link clinical variables with genomic alterations.
- Locate exact study IDs and clinical field names directly in chat without clicking through web filters.

### Research & Workflow Use Cases

- Estimate mutation frequencies and hotspot distributions for candidate genes across specific tumor types.
- Compare sample composition and clinical covariates across studies before pooling cohorts for downstream survival or association analyses.
- Cross-check genomic alteration patterns and clinical annotations when prioritizing therapeutic targets or biomarkers.
- Verify cohort metadata, sample counts, and primary study citations during manuscript or grant preparation.

### Usage Notes

- Sign in via OAuth on first connection. This official cBioPortal MCP endpoint queries published datasets and does not bypass controlled-access data agreements.
- Available studies, clinical variables, and genome builds reflect the upstream cBioPortal release at the time of query.
- Homepage: [cBioPortal MCP Documentation](https://docs.cbioportal.org/ai-integrations/mcp/). Cite both cBioPortal and the originating study publications when reporting results.
