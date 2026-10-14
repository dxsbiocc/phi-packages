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
| HGNC | complete set | Approved symbols, previous symbols, aliases, Entrez ids, and Ensembl ids |
| Cellosaurus | release 56.0 | Accession, name, synonym, and problematic-cell-line comments |

DepMap's catalogue endpoints `https://depmap.org/portal/api/download/files` and `https://depmap.org/portal/api/no-captcha/download/files` returned an HTTP 200 Cloudflare HTML page from the network used to build this pipeline. Figshare article records answered, so `sources.yaml` stores those download URLs. Each DepMap and PRISM source also lists resolvers. `phi-pgx download --resolve-urls` tries the portal catalogue first and, if that body is HTML or has no URL, uses the Figshare article. A file's `url` is the reproducible default.

The downloader writes `raw/<source>/<version>/` and `sha256sum.json`. A finished file whose hash matches the manifest is skipped. A `.partial` file is resumed with HTTP `Range` when the server honors it. An HTML body is rejected, so a bot-check page is not stored as a CSV. Optional `sha256` values in `sources.yaml` fail the download on mismatch.

## Identifiers

Normalization in `pgx.normalize` runs before identifier mapping and is shared by DepMap, PRISM, GDSC, and Cell Model Passports. One module handles each entity. Unmatched and ambiguous rows are written to `review/` and are not guessed.

### Genes

DepMap headers such as `TP53 (7157)` are split into a source symbol and an Entrez id. A header can also carry no id, several ids (`7157;7158` or `7157, ENSG00000141510`), or a token that is not a gene symbol. The HGNC complete set (`sources.yaml`, `hgnc`) maps an Entrez id, then an approved symbol, then a previous symbol, then an alias. The approved symbol is stored in `symbol`. The header symbol stays in `source_symbol`, with `match_method` set to `hgnc_entrez`, `hgnc_symbol`, `hgnc_previous`, or `hgnc_alias`. An HGNC Ensembl id is kept when the header does not disagree with it. Several Entrez ids, or a symbol and an Entrez id that point at different genes, are `ambiguous`.

### Cell lines

Cell lines are keyed by a Cellosaurus accession (`CVCL_` plus four letters or digits). Names are case-folded and stripped of punctuation, spaces, hyphens, and underscores, so `NCI-H460` and `NCIH460` are the same key. `H460` is not: the normalizer does not drop an `NCI` prefix. `H460` matches NCI-H460 only when Cellosaurus lists it as a synonym (release 56 does, as `CVCL_0459`). A DepMap `CCLEName` tissue suffix such as `_LUNG` is removed for a second lookup key; the original name is kept. An explicit RRID wins over a name. A name or synonym is used only when the row has no RRID, SIDM, COSMIC, or ACH id that failed to resolve, and only when that key hits exactly one accession. Cellosaurus `Problematic cell line` comments set `cellosaurus_flag` to `misidentified`, `contaminated`, or `problematic` and the row stays mapped. Two accessions, or an RRID that disagrees with the name, stay `ambiguous`.

Without `cellosaurus.txt`, name matching falls back to Cell Model Passports names, with the same rule that a stable id blocks the name fallback.

### Compounds

Compounds are keyed by PubChem CID. Source names are whitespace- and Unicode-normalized for lookup. Stereo and enantiomer parentheticals (`Nutlin-3a (-)`) and trailing salt tokens are removed only from the lookup string; `name` keeps the source text. Synonym lists are split on `;`, `|`, and commas. GDSC `DRUG_ID` is stored as `gdsc_drug_id` and is not sent to PubChem as a name. BRD ids are a last lookup when the name and any SMILES or InChIKey missed.

`phi-pgx resolve-compounds` prefers a source InChIKey, then SMILES, then the exact name, then the stereo-stripped name, then the salt-stripped name. The first query that returns one CID wins and records `match_method` plus `confidence` (`high` for structure, `medium` for an exact name or BRD, `low` for a stripped name). More than one CID is `ambiguous` and a looser query is not tried. ChEMBL ids come from a unique standard InChIKey match. Lookups are cached under `raw/pubchem/`, rate-limited, and resumed from that cache. A multi-name PUG body is not used, because it returns one CID list and drops which name matched.

## Tables

Processed output:

- `model_xref/` and `compound_xref/`
- `gene_effect/`, `expression/`, `cnv/`, `mutation/`, `drug_response/` partitioned as `dataset=<name>/`
- `releases/` with version, license, and citation
- `gene_xref/` when the HGNC complete set and a gene matrix or mutation file are present
- `review/models/`, `review/compounds/`, and `review/genes/` for unresolved primary ids
- `validation/report.json` after `phi-pgx validate`

`drug_response.metric` is the published statistic: GDSC `AUC` and `LN_IC50` (natural log of micromolar; lower means a stronger response), PRISM secondary `AUC`, `IC50`, and `EC50`, and PRISM primary / 24Q2 `LFC`. Expression is `log2_tpm_plus_1` from the DepMap `TPMLogp1` matrix, not raw TPM. Copy number is the relative gene-level value in `OmicsCNGene.csv`. Mutations are the DepMap hg38 calls.

`validate` checks required columns and writes mapping coverage. The Nutlin-3a check expects a lower median GDSC `LN_IC50` and `AUC` in TP53 wild-type lines than in lines with a non-silent `ncbigene:7157` mutation. Groups smaller than 3 models are skipped. Lines absent from the mutation table are excluded. A failed check exits non-zero.

## Later batches

`sources.yaml` `planned` lists sources that are not downloaded or parsed:

- NCI60 DTP April 2026
- Frozen ORCESTRA PSets CTRPv2.1, gCSI, FIMM, GRAY, UHNBreast, and CCLE, via R/PharmacoGx or rpy2
- HMDB bulk XML

Do not point this pipeline at an ORCESTRA PRISM PSet. PRISM stays on the DepMap Figshare releases above.
