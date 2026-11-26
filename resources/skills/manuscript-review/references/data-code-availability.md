# FAIR Data & Code Availability Branch

Load this reference ONLY when planning data/code deposition, auditing human-data privacy boundaries, or drafting Data & Code Availability statements.

## 1. Domain Repository Routing Matrix

| Data Type | Open-Access Repositories | Controlled-Access (Human Identifiable) |
| --- | --- | --- |
| **Bulk / scRNA-seq / Spatial / Epigenomics** | **GEO** / **ArrayExpress** (processed) + **SRA** / **ENA** / **GSA** (raw) | **EGA**, **dbGaP**, **GSA-Human** (raw reads + sensitive clinical metadata) |
| **WGS / WES / Variants** | **ENA** / **SRA** / **EVA** (non-human / consented open) | **EGA**, **dbGaP**, **GSA-Human** |
| **Proteomics / Metabolomics** | **PRIDE** / **MassIVE** / **iProX**; **MetaboLights** / **Metabolomics Workbench** | Institutional DAC + controlled portal |
| **Microscopy / Bioimaging** | **BioImage Archive (BIA)**, **IDR**, **EMPIAR** | De-identified DICOM / controlled archive |
| **Tabular Source Data / Clinical** | **Zenodo**, **Dryad**, **Figshare**, Journal Source Data | **Vivli**, **YODA**, **dbGaP** |
| **Analysis Code & Notebooks** | **GitHub/GitLab** + **versioned Zenodo / Software Heritage DOI** | Same (never commit PHI or API keys!) |

## 2. Hard Compliance Rules

- **Privacy Hard Stop**: Never deposit identifiable human raw reads (`FASTQ/BAM/VCF`), PHI, or regulated human genetic resources into open repositories; require controlled-access repository + IRB ID + DAC contact & response timeline.
- **No Generic "Upon Reasonable Request"**: Flag as non-compliant unless backed by documented legal/ethical restriction + DAC contact.
- **Zero Fake Accessions**: Use `[GEO GSE###### TBD]`, `[HRA###### TBD]`, `[ZENODO DOI TBD]` until minted.
- **Workspace Reproducibility Scan**: Check for `environment.yml`/`renv.lock`/`Dockerfile`, `LICENSE`, random seeds, and scrubbed secrets/PHI.
