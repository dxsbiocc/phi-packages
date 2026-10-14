# Pharmacogenomics pipeline

Batch 1 turns DepMap, PRISM, GDSC, and Cell Model Passports releases into partitioned Parquet. The MCP server in `mcp-servers/pharmacogenomics` reads that directory. It does not download the releases itself.

```bash
pip install -e '.[test]'
pytest
phi-pgx download          # small metadata and crosswalks
phi-pgx download --live   # CRISPR, expression, mutations, copy number, dose-response matrices
phi-pgx process
phi-pgx validate
PHI_PGX_DATA_DIR=processed phi-omics-pharmacogenomics
```

`download` without `--live` skips files marked `live: true`. Those matrices are hundreds of megabytes to about 800 MB. `PHI_PGX_LIVE=1` is the same switch. Tests use fixtures and do not download them.

## Sources

| Source | Version | What batch 1 uses |
| --- | --- | --- |
| DepMap | 24Q2 | `Model.csv`, `CRISPRGeneEffect.csv`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, `OmicsCNGene.csv`, `OmicsSomaticMutations.csv` |
| PRISM | 19Q4 | Primary-screen LFC and secondary-screen dose-response curve parameters, plus cell and compound metadata. This is the Corsello 2020 release on Figshare, not an ORCESTRA PRISM PSet. |
| PRISM | Public 24Q2 | Collapsed LFC and compound metadata from the DepMap Figshare record |
| GDSC | release 8.4 | `screened_compounds_rel_8.4.csv`, GDSC1 and GDSC2 fitted dose response dated 24 Jul 2022. URLs are pinned to `release-8.4`, which is the same tree `current_release` pointed at when this file was written. |
| Cell Model Passports | model list | `model_list_latest.csv` for SIDM, BROAD/ACH, COSMIC, and RRID |

DepMap's catalogue endpoints `https://depmap.org/portal/api/download/files` and `https://depmap.org/portal/api/no-captcha/download/files` returned an HTTP 200 Cloudflare HTML page from the network used to build this pipeline. Figshare article records answered, so `sources.yaml` stores those download URLs. Each DepMap and PRISM source also lists resolvers. `phi-pgx download --resolve-urls` tries the portal catalogue first and, if that body is HTML or has no URL, uses the Figshare article. A file's `url` is the reproducible default.

The downloader writes `raw/<source>/<version>/` and `sha256sum.json`. A finished file whose hash matches the manifest is skipped. A `.partial` file is resumed with HTTP `Range` when the server honors it. An HTML body is rejected, so a bot-check page is not stored as a CSV. Optional `sha256` values in `sources.yaml` fail the download on mismatch.

## Identifiers

Cell lines are keyed by a Cellosaurus accession (`CVCL_` plus four letters or digits). A field that lists two accessions is ambiguous. DepMap `RRID`, Cell Model Passports `RRID`, Sanger `SIDM`, and COSMIC ids are explicit crosswalks. A normalized model name is used only when the row has no RRID, SIDM, or COSMIC id and that name matches exactly one accession. If two accessions disagree, or a name matches several, the row is `ambiguous` and `cvcl` is null. Nothing is fuzzy-matched.

Compounds are keyed by PubChem CID, with ChEMBL id and InChIKey when the source or a unique PubChem name lookup provides exactly one. GDSC 8.4 `screened_compounds_rel_8.4.csv` has names and drug ids, not CIDs. PRISM files have BRD ids, names, and sometimes SMILES. Those rows stay `unmatched` until `phi-pgx resolve-compounds`, which calls PubChem for a unique name and records the method. Ambiguous name hits stay unresolved. Source names, BRD ids, and GDSC drug ids are kept on `compound_xref`.

Genes are `ncbigene:<entrez>` plus the HGNC symbol parsed from labels such as `TP53 (7157)`.

## Tables

Processed output:

- `model_xref/` and `compound_xref/`
- `gene_effect/`, `expression/`, `cnv/`, `mutation/`, `drug_response/` partitioned as `dataset=<name>/`
- `releases/` with version, license, and citation
- `review/models/` and `review/compounds/` for unresolved primary ids
- `validation/report.json` after `phi-pgx validate`

`drug_response.metric` is the published statistic: GDSC `AUC` and `LN_IC50` (natural log of micromolar; lower means a stronger response), PRISM secondary `AUC`, `IC50`, and `EC50`, and PRISM primary / 24Q2 `LFC`. Expression is `log2_tpm_plus_1` from the DepMap `TPMLogp1` matrix, not raw TPM. Copy number is the relative gene-level value in `OmicsCNGene.csv`. Mutations are the DepMap hg38 calls.

`validate` checks required columns and writes mapping coverage. The Nutlin-3a check expects a lower median GDSC `LN_IC50` and `AUC` in TP53 wild-type lines than in lines with a non-silent `ncbigene:7157` mutation. Groups smaller than 3 models are skipped. Lines absent from the mutation table are excluded. A failed check exits non-zero.

## Later batches

`sources.yaml` `planned` lists sources that are not downloaded or parsed:

- NCI60 DTP April 2026
- Frozen ORCESTRA PSets CTRPv2.1, gCSI, FIMM, GRAY, UHNBreast, and CCLE, via R/PharmacoGx or rpy2
- HMDB bulk XML

Do not point this pipeline at an ORCESTRA PRISM PSet. PRISM stays on the DepMap Figshare releases above.
