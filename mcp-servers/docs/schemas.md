# Shared omics schemas

Every tool on every server returns one JSON object. It is either a page or an error. Source-specific fields that do not fit the shared entity go in `source_extra` and nowhere else.

The Zod schemas and inferred TypeScript types live in `mcp-servers/schemas`. Import them from the workspace package `phi-omics-schemas` (the package is private until an npm scope is chosen).

## Identifiers

Identifiers are CURIEs: `prefix:local`. Dotted prefixes are allowed (`gdc.case`, `gtex.tissue`, `cellxgene.dataset`).

| Concept | CURIE |
| --- | --- |
| Ensembl gene | `ensembl:ENSG00000141510` (version stripped; the versioned GENCODE id stays in `source_extra.gencode_id`) |
| NCBI gene | `ncbigene:7157` |
| Symbol only | `hgnc.symbol:TP53` when no stable database id is available |
| UniProt | `uniprot:P04637` |
| Taxon | `NCBITaxon:9606` |
| Uberon, CL, EFO, MONDO, PATO, CHEBI | ontology prefix plus zero-padded local id, for example `UBERON:0002190` or `CL:0000636` |
| Source records | `gtex.tissue:Whole_Blood`, `gdc.case:<uuid>`, `gdc.file:<uuid>`, `gdc.ssm:<uuid>`, `encode:ENCSR...`, `pride:PXD001357`, `pride.file:<id>`, `mtbls:MTBLS1`, `hubmap:HBM...`, `hca.project:<uuid>`, `cellxgene.collection:<uuid>`, `cellxgene.dataset:<uuid>`, `cbioportal.study:<id>`, `pdc.study:PDC000759` |
| Cell line | `cellosaurus:CVCL_0465` |
| Compound | `pubchem:11433190`, with `chembl:CHEMBL191334` in `xrefs` when a single ChEMBL id is known |
| Drug combination | `gdsc.combination:<anchor>:<library>` or `gdsc.combination:annotation:<drug id>`. This is not a compound id |
| Drug response | `pgx.response:<dataset>/<model>/<compound>/<metric>` |
| Pharmacogenomics dataset | `pgx.dataset:gdsc1` |
| Publication | `doi:<lowercased bare DOI>` when a DOI exists, even if there is no PMID. Otherwise `pmid:<digits>`, then `arxiv:<id without version>`. `pmcid:PMC…`, `openalex:W…`, and `s2:<id>` are xrefs, not primary keys |

`xrefs` is an array of additional CURIEs on every entity. It is always present.

## Coordinates and expression

`locus` is null or:

- `genome_build`: string such as `GRCh38` or `hg19`
- `chrom`: `chr`-prefixed (`chr17`, `chrX`, `chrM`)
- `start`, `end`: 1-based inclusive integers
- `strand`: `+`, `-`, `.`, or null

Expression `unit` is one of `TPM`, `log2FC`, `count`, `UMI`, `protein_level`, `z_score`. `statistic` names the summary (`median`, `log2_allelic_fold_change`). GTEx eGene fold changes use `log2FC`. GTEx median expression uses `TPM`.

## Entities

Each entity has `entity_type`, `id`, `xrefs`, and `source_extra`, plus the fields below. Missing values are `null` or empty arrays, not omitted keys.

| entity_type | Fields |
| --- | --- |
| `gene` | `symbol`, `name`, `taxon_id`, `biotype`, `synonyms`, `locus` |
| `sample` | `name`, `donor_id`, `tissue_id`, `cell_type_id`, `disease_ids`, `assay`, `taxon_id` |
| `donor` | `sex` (`female`, `male`, `unknown`, or null), `age_years`, `taxon_id`, `disease_ids` |
| `dataset` | `title`, `description`, `assay`, `taxon_id`, `tissue_ids`, `disease_ids`, `cell_type_ids`, `sample_count`, `access`, `url` |
| `file` | `name`, `format`, `data_type`, `size_bytes`, `access`, `url`, `dataset_id` |
| `cell_type` | `name`, `description`, `ontology_id`, `parent_ids` |
| `tissue` | `name`, `ontology_id` |
| `protein` | `name`, `gene_id`, `gene_symbol`, `taxon_id` |
| `metabolite` | `name`, `formula`, `inchikey` |
| `compound` | `name`, `pubchem_cid`, `chembl_id`, `inchikey`, `synonyms` |
| `drug_combination` | `name`, `components` (`compound_id`, `role` `anchor` or `library`, `name`, `concentrations`, `unit`) |
| `drug_response` | `compound_id`, `combination_id`, `sample_id`, `dataset_id`, `metric`, `value`, `unit` |
| `expression` | `gene_id`, `tissue_id`, `cell_type_id`, `sample_id`, `value`, `unit`, `statistic` |
| `variant` | `gene_id`, `locus`, `consequence` |
| `publication` | `doi`, `title`, `authors` (`name`, `orcid`), `venue`, `dates` (`published`, `posted`, `updated`), `abstract`, `type`, `is_oa`, `oa_urls`, `license`, `published_as`, `preprint_of` |

`publication.type` is `journal-article`, `preprint`, `review`, `book-chapter`, `proceedings-article`, `other`, or null. `oa_urls[]` has `url`, `content_type` (`html`, `pdf`, `xml`, `landing`, `other`), and `version`. A preprint and its journal version are two records. `published_as` on the preprint and `preprint_of` on the journal article point at the other record's CURIE. When `doi` is set it is the lowercased DOI without `https://doi.org/`, and `id` is `doi:` plus that value. PMID and PMCID then live in `xrefs`.

`access` is `open`, `controlled`, or `unknown`. A controlled file must have `url: null`. The schema rejects any other value. Controlled datasets also omit portal URLs when the URL would point at the data rather than an open landing page; this codebase sets `url` to null whenever `access` is `controlled`.

`metabolite` is part of the shared model. None of the verified MetaboLights JSON endpoints return a parsed metabolite table, so the current tools do not emit metabolite rows.

`compound`, `drug_combination`, and `drug_response` are the pharmacogenomics entities. A compound id is `pubchem:<CID>` when one PubChem CID is known. A combination is a different entity: its id is `gdsc.combination:<anchor>:<library>`, and `components` lists each drug in anchor then library order with that drug's own compound id. A combination name is never stored as a compound. `drug_response.combination_id` is null for a single agent and set for a combination result, whose `compound_id` stays null. `metric` is the published statistic (`AUC`, `IC50`, `LN_IC50`, `EC50`, `LFC`, or a combination delta such as `SYNERGY_DELTA_AUC`). `unit` is null when the release does not name one. DepMap expression is log2(TPM+1), so those rows use `unit: null` and `statistic: "log2_tpm_plus_1"` rather than `TPM`.

## Page

```json
{
  "entity_type": "page",
  "item_type": "gene",
  "items": [],
  "page": { "offset": 0, "limit": 20, "returned": 0, "total": 0, "next_offset": null },
  "source": { "id": "gtex", "name": "GTEx Portal", "version": "gtex_v10", "homepage": "https://gtexportal.org/", "retrieved_at": "2026-10-14T00:00:00.000Z" },
  "usage": { "summary": "...", "url": "https://gtexportal.org/home/datasets" },
  "source_extra": {}
}
```

`item_type` matches every `items[].entity_type`. `page.returned` is `items.length`. `next_offset` is null on the last page. Tool arguments use the same `offset` and `limit` (default 20, maximum 50).

`usage_terms` returns an empty page. The text is in `usage.summary`.

## Error

```json
{
  "entity_type": "error",
  "code": "upstream_http",
  "message": "...",
  "source": { "id": "gdc", "name": "NCI GDC", "version": "GRCh38", "homepage": "https://gdc.cancer.gov/", "retrieved_at": "2026-10-14T00:00:00.000Z" },
  "retryable": true,
  "source_extra": { "status": 503 }
}
```

`code` is `upstream_http`, `not_found`, `invalid_input`, or `unavailable`.
