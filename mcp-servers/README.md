# Omics MCP servers

Seven read-only MCP servers for public scientific databases that do not ship an official MCP endpoint. They share one entity model (`docs/schemas.md`). The remote servers share one HTTP client (timeouts, retries, per-host rate limit, on-disk cache). Pharmacogenomics reads local Parquet instead.

| Server | Package | Sources |
| --- | --- | --- |
| Transcriptomics and single-cell | `phi-omics-transcriptomics` | GTEx v10, CELLxGENE Discover, CellGuide, HCA Azul, HuBMAP search |
| Cancer omics | `phi-omics-cancer` | GDC, cBioPortal public REST, PDC GraphQL |
| Epigenomics | `phi-omics-epigenomics` | ENCODE portal, Roadmap via `award.rfa=Roadmap` |
| Proteomics | `phi-omics-proteomics` | PRIDE, Human Protein Atlas |
| Metabolomics | `phi-omics-metabolomics` | MetaboLights |
| Pharmacogenomics | `phi-omics-pharmacogenomics` | DepMap 24Q2, PRISM 19Q4 and Public 24Q2, GDSC 8.4, Cell Model Passports |
| Literature | `phi-literature` | PubMed, Europe PMC, bioRxiv/medRxiv, Crossref, OpenAlex, Semantic Scholar, arXiv, Unpaywall, OpenCitations |

The same TypeScript runs locally over stdio or Streamable HTTP, and in Docker on one host. Packages are private workspace members. An npm scope and a publish workflow are intentionally not set up yet; `bin`, `files`, and `engines` are already in each package so that step can be added later without moving code.

Public HTTPS endpoints for the five omics servers are live on `phi.omicsagent.cn`. Matching manifests live under `resources/connectors/` (`transcriptomics-single-cell`, `cancer-omics`, `epigenomics`, `proteomics`, `metabolomics`). Pharmacogenomics is mounted at `/pharmacogenomics/mcp` and reads local Parquet. Literature is mounted at `/literature/mcp`. Add those paths to the catalog only after they are deployed.

## Build and test

From the repository root (npm workspaces):

```bash
npm install
npm run build:mcp
npm run lint:mcp
npm test:mcp
```

`npm run test:mcp` runs the TypeScript tests, including live API checks. Set `PHI_MCP_SKIP_LIVE=1` to skip those checks. Content catalog tests stay separate:

```bash
node --test tests/skills.test.mjs tests/content-boundaries.test.mjs tests/content-icons.test.mjs
```

## Local use

Each server defaults to stdio, which is what desktop MCP clients spawn. HTTP mode binds to `127.0.0.1` unless `HOST` is set.

```bash
node mcp-servers/transcriptomics-single-cell/dist/cli.js
node mcp-servers/transcriptomics-single-cell/dist/cli.js --http --port 8787
node mcp-servers/cancer-omics/dist/cli.js --help
npm exec -w phi-omics-epigenomics -- phi-omics-epigenomics --stdio
npm exec -w phi-omics-proteomics -- phi-omics-proteomics --http --port 8788
npm exec -w phi-omics-metabolomics -- phi-omics-metabolomics --http --port 8789
PHI_PGX_DATA_DIR=/path/to/processed npm exec -w phi-omics-pharmacogenomics -- phi-omics-pharmacogenomics --http --port 8790
npm exec -w phi-literature -- phi-literature --http --port 8791
```

The bins are `phi-omics-transcriptomics`, `phi-omics-cancer`, `phi-omics-epigenomics`, `phi-omics-proteomics`, `phi-omics-metabolomics`, `phi-omics-pharmacogenomics`, `phi-literature`, and `phi-omics-gateway`. Each server defaults to stdio and switches to HTTP with `--http`. The remote servers and literature have no native addons. Pharmacogenomics uses the DuckDB Node binding (`@duckdb/node-api`) to read Parquet.

Example stdio entry for a desktop client, with paths adjusted to this checkout:

```json
{
  "mcpServers": {
    "transcriptomics": {
      "command": "node",
      "args": ["/absolute/path/phi-packages/mcp-servers/transcriptomics-single-cell/dist/cli.js"]
    }
  }
}
```

HTTP endpoints for one server are `http://127.0.0.1:8787/mcp`, `GET /health`, and `GET /ready`.

Cache files go to `PHI_OMICS_CACHE_DIR`, or `$XDG_CACHE_HOME/phi-omics-mcp`, or `~/.cache/phi-omics-mcp`. On Windows the default is `%LOCALAPPDATA%\phi-omics-mcp\cache`.

| Variable | Default | Purpose |
| --- | --- | --- |
| `HOST` | `127.0.0.1` | HTTP bind address. Use `0.0.0.0` in Docker. |
| `PORT` | `8080` | HTTP port |
| `BASE_PATH` | empty | Prefix in front of `/health` and `/mcp` |
| `PHI_OMICS_CACHE_DIR` | platform cache dir | GET/POST response cache |
| `PHI_OMICS_CACHE_TTL_MS` | 6 hours | Cache lifetime |
| `PHI_OMICS_RATE_LIMIT_RPS` | `4` | Per-host request spacing |
| `PHI_OMICS_TIMEOUT_MS` | `30000` | Upstream timeout |
| `PHI_OMICS_RETRIES` | `2` | Retries for network errors, 429, and 5xx |
| `PHI_OMICS_MAX_BODY_BYTES` | `8000000` | Reject larger upstream bodies while they stream, before they are fully buffered |
| `PHI_OMICS_CACHE_MAX_BYTES` | `500000000` | Cache size cap. Expired files are removed first, then the oldest, at most every ten minutes |
| `PHI_OMICS_INBOUND_RPM` | `0` (off) | HTTP mode: requests per minute allowed from one client on the MCP paths. Over the limit returns 429 with `Retry-After` |
| `PHI_OMICS_TRUST_PROXY` | unset | Set to `1` behind a reverse proxy so the inbound limit keys on the first `X-Forwarded-For` address. Leave unset when clients can reach the process directly |
| `PHI_OMICS_CORS_ORIGIN` | `*` | `Access-Control-Allow-Origin` value. Set the site origin when browsers are the only clients |
| `PHI_OMICS_ACCESS_LOG` | on | Set to `0` to stop the one-line JSON request log written to stderr. Tool errors are always logged |
| `GTEX_DATASET` | `gtex_v10` | GTEx dataset |
| `GTEX_GENCODE` | `v39` | GTEx GENCODE version |
| `GTEX_GENOME` | `GRCh38/hg38` | GTEx genome build argument |
| `CELLGUIDE_SNAPSHOT` | `1692975365` | Pinned CellGuide public snapshot |
| `MCP_SERVER` | `all` | Gateway only: one server name, or `all` |
| `PHI_PGX_DATA_DIR` | unset | Pharmacogenomics only: directory of processed Parquet from `data-pipelines/pharmacogenomics` |
| `NCBI_API_KEY` | empty | Literature: NCBI E-utilities key. Without it, keep the shared rate near 3 requests/second |
| `PHI_LITERATURE_EMAIL` | empty | Literature: contact email for the Crossref, OpenAlex, and Unpaywall polite pools. Also read from `NCBI_EMAIL` or `UNPAYWALL_EMAIL` |
| `S2_API_KEY` | empty | Literature: Semantic Scholar key, used only when that source is requested |

No upstream API key is required. The servers do not send tokens and do not download controlled files.

## One host with Docker

From `mcp-servers/`:

```bash
docker compose up --build
```

The image listens on port 8080:

- `GET /health` (liveness: the process is up)
- `GET /ready` (readiness: the cache directory is writable; each server also reports whether it has what it needs, such as the pharmacogenomics data directory. A server missing its data is listed with a reason but does not fail the probe)
- `POST /transcriptomics/mcp`
- `POST /cancer/mcp`
- `POST /epigenomics/mcp`
- `POST /proteomics/mcp`
- `POST /metabolomics/mcp`
- `POST /pharmacogenomics/mcp`
- `POST /literature/mcp`

`MCP_SERVER=transcriptomics-single-cell` (or `cancer-omics`, `epigenomics`, `proteomics`, `metabolomics`, `pharmacogenomics`, `literature`) serves that one server at `/mcp` instead.

Mount the processed pharmacogenomics tables by setting `PHI_PGX_DATA_DIR` on the host to that directory. Compose maps it read-only to `/var/lib/phi-pgx`. The pipeline itself runs outside the image.

Put TLS on a reverse proxy. `Caddyfile.example` proxies a hostname to `127.0.0.1:8080`. Caddy obtains and renews certificates. Nginx or another proxy is the same pattern: forward `POST`, `GET`, and `DELETE` to the Node process, preserve `Authorization` and `Mcp-Session-Id` if you add auth later, and do not buffer streaming responses. This process does not terminate TLS and does not authenticate clients. It can limit request rates per client (`PHI_OMICS_INBOUND_RPM`), but a public deployment should also rate-limit at the proxy. Bind it to localhost when the proxy is on the same machine.

## Tools

Pagination arguments on list tools are `offset` (default 0) and `limit` (default 20, maximum 50). Every tool returns the shared page or error object.

### Transcriptomics and single-cell

`usage_terms`, `gtex_search_genes`, `gtex_list_tissues`, `gtex_gene_expression`, `gtex_list_egenes`, `gtex_list_eqtls`, `cellxgene_list_collections`, `cellxgene_get_collection`, `cellxgene_list_datasets`, `cellxgene_get_dataset`, `cellguide_cell_type`, `hca_search_projects`, `hca_get_project`, `hubmap_search_datasets`, `hubmap_search_samples`, `hubmap_search_donors`.

CELLxGENE Census matrices are not exposed. The Census reader is a native TileDB library, and this server stays on the Discover curation JSON API. `cellxgene_get_dataset` reads the dataset out of its parent collection because `GET /curation/v1/datasets/{id}` returned 404.

### Cancer omics

`usage_terms`, `gdc_search_cases`, `gdc_get_case`, `gdc_search_files`, `gdc_get_gene`, `gdc_search_mutations`, `cbioportal_list_studies`, `cbioportal_get_study`, `cbioportal_list_samples`, `cbioportal_list_patients`, `cbioportal_get_gene`, `pdc_list_studies`, `pdc_get_study`.

cBioPortal here is the public REST API. The existing `cbioportal` connector (OAuth, `https://mcp.cbioportal.org/db/mcp`) is a different service and is left in place.

### Epigenomics

`usage_terms`, `encode_search_experiments`, `encode_get_experiment`, `encode_search_roadmap`, `encode_search_files`, `encode_get_file`, `encode_search_biosamples`, `encode_get_biosample`.

ENCODE's edge blocks the `from` query parameter, so search tools request the leading window and slice it. `offset + limit` must be at most 100.

### Proteomics

`usage_terms`, `pride_search_projects`, `pride_get_project`, `pride_list_files`, `hpa_search_genes`, `hpa_get_gene`, `hpa_list_proteins`.

### Metabolomics

`usage_terms`, `metabolights_list_studies`, `metabolights_get_study`, `metabolights_list_files`, `metabolights_list_assays`.

### Pharmacogenomics

`usage_terms`, `resolve_cell_line`, `list_tissues`, `gene_dependency`, `gene_expression`, `gene_copy_number`, `gene_mutations`, `compound_response`, `list_response_metrics`, `biomarker_association`, `model_summary`.

- **Tissue filter.** `gene_dependency`, `gene_expression`, `gene_copy_number`, `gene_mutations`, `compound_response` and `biomarker_association` accept `tissue`, matched case-insensitively against the tissue of mapped Cellosaurus models. `list_tissues` returns the values. An unknown tissue is a `not_found` error, never an unfiltered result.
- **Gene tools.** `gene_expression` and `gene_copy_number` order rows by value (`order=asc|desc`) and put the distribution (n, median, quartiles, range) of every matching row in `source_extra.distribution`. `gene_mutations` leaves silent variants out unless `include_silent` is true, filters by `consequence` and `protein_change`, and reports `mutated_models`, `sequenced_models`, `frequency` and the top protein changes. A model counts as wild type only if it appears in the mutation table.
- **Response metrics.** The same metric name does not mean the same thing in every dataset: Sanger GDSC and PRISM publish an AUC that falls as a model gets more sensitive, while an AUC recomputed by PharmacoGx for an ORCESTRA PharmacoSet may be a response area that rises with sensitivity. `list_response_metrics` lists every dataset and metric with row, model and compound counts and a `direction` (`lower_is_stronger`, `higher_is_stronger`, or `unknown`). Direction is only asserted where the source documents it (potency metrics such as `LN_IC50`, `IC50`, `GI50`, `LFC` and `VIABILITY`; `AAC`; `AUC` for `gdsc1`, `gdsc2`, `prism_secondary_19q4` and `ctrpv2`). It is `unknown` for PharmacoSet `AUC` and combination deltas rather than guessed. `compound_response` groups rows by dataset and metric and orders by value only within a group, because values from different groups are not comparable. Each row carries `direction`, the page lists `metric_groups` (rows, models, median, range), and the summary warns when rows mix groups. Set `dataset` and `metric` to compare models. Earlier versions ordered all rows by value, which interleaved incomparable metrics. When the processed data was built with the standardization step, each row also carries `sensitivity_percentile` (among all rows of the same dataset, metric and dose) and `compound_percentile` (among rows of the same compound), both from 0 least sensitive to 1 most sensitive and null where the direction is `unknown`; older parquet returns null for both. IC50, EC50 and AC50 can also carry `beyond_tested_range`: true when the value is larger than the highest dose tested for that model and compound, so it was extrapolated from an incomplete curve and should be treated as censored rather than as a measurement (the value is returned unchanged; AUC or AAC do not have this problem). It is null when no raw doses were loaded.
- **`biomarker_association`.** For one single-agent compound and one gene, it tests mutation status (Mann-Whitney U, both groups at least `min_group` models, default 3) or expression, copy number, or CRISPR dependency (Spearman, at least 10 paired models), once per response dataset, metric and dose. Rows are sorted by p-value and carry medians or rho, n, and a rank-biserial effect size. When the metric's direction is known, a row also says which group (`more_sensitive_group`) or feature level (`higher_feature_more_sensitive`) responds more strongly; it is null when the direction is unknown. Groups that are too small are listed in `source_extra.skipped_groups` with the reason. The statistics are descriptive: p-values are uncorrected for multiple testing, cell lines are not independent patients, and tissue confounds many hits, so repeat with `tissue` set.

`compound_response` returns single-agent rows unless `combinations` is `include` or `only`, or the query is a combination id or an `Anchor + Library` name. Combination rows set `combination_id` and leave `compound_id` null.

The server does not download releases. Build the tables with `data-pipelines/pharmacogenomics` and point `PHI_PGX_DATA_DIR` at the `processed/` directory. Without that directory the query tools return `unavailable`. `usage_terms` still returns the license text.

### Literature

`usage_terms`, `search_publications`, `get_publication`, `get_fulltext`, `related_publications`, `resolve_references`.

The tools follow the jobs a person actually does with the literature, not one tool per website. Search hits are short: at most eight authors and a 400-character abstract unless `include_abstract` is true. `get_publication` returns the abstract. Full text is inlined only when a Creative Commons or public-domain license is known; otherwise the response is links.

| Tool | Task | Why this shape |
| --- | --- | --- |
| `usage_terms` | Same contract as the other servers | Names the routing rules and the optional email and API key |
| `search_publications` | Find papers on a topic, gene, or disease, newest first. Also the newest preprints on a topic, or a date window of bioRxiv/medRxiv | PubMed is the default because MeSH and field tags are how biomedical search stays precise. `type=preprint` uses Europe PMC. An empty query plus `server` and dates browses bioRxiv or medRxiv, which is what that API actually supports |
| `get_publication` | Open one paper from a DOI, PMID, PMCID, arXiv id, or a title or citation string. See whether a preprint was published | One record carries the abstract, ORCID authors, license, open links, and `published_as` / `preprint_of`. A separate id-converter or copyright call is not needed |
| `get_fulltext` | Read or quote the methods from the open copy | Europe PMC JATS is returned as clipped sections when the license allows it. A DOI that OpenAlex does not tag with a PMCID is looked up in Europe PMC before that fetch. PDFs are links, not parsed text |
| `related_publications` | Move from a key paper to similar papers, citations, or references | `relation` picks one of those. PubMed supplies similar articles when a PMID exists. OpenAlex supplies citations, references, and DOI-only similarity. OpenCitations is the fallback, not a second call in the happy path |
| `resolve_references` | Turn a reference list into records | Up to 25 strings. DOIs and PMIDs are batched. A line that does not resolve is listed in `source_extra.unresolved` |

Routing is deliberate. A biomedical query goes to PubMed. A preprint topic goes to Europe PMC. A DOI goes to OpenAlex, then Crossref if OpenAlex has no work. A `10.1101` DOI goes to the bioRxiv/medRxiv details API so the journal DOI can be returned in `published_as`. arXiv and Semantic Scholar run only when `sources` or the id asks for them. Passing `sources` queries at most three of those databases and merges duplicates on DOI, then PMID, then PMCID.

Identifiers: `doi:` plus the lowercased DOI with any `https://doi.org/` prefix removed, including online-first articles that have no PMID. If there is no DOI, the id is `pmid:`, then `arxiv:`. PMCID, OpenAlex, and Semantic Scholar ids stay in `xrefs`. The preprint and the journal version are two records.

#### How this compares with other literature connectors

Designs were read from public docs and READMEs. No code was copied.

| Reference | Their tools | What they get right | Gap this server does not repeat |
| --- | --- | --- | --- |
| Anthropic PubMed connector (`pubmed.mcp.claude.com`, seven tools: `search_articles`, `get_article_metadata`, `find_related_articles`, `lookup_article_by_citation`, `convert_article_ids`, `get_full_text_article`, `get_copyright_status`) | Search, one record, related articles, citation lookup, id conversion, PMC full text, copyright | Task-shaped and small. PubMed query syntax. Citation lookup is its own step. Full text is separate from the abstract | PubMed and PMC only, so a DOI-only online-first article and a preprint are out of scope. Id conversion and copyright are extra hops. Related articles are not the citing or cited lists |
| Anthropic bioRxiv/medRxiv connector (nine tools, including `search_preprints`, `get_preprint`, `search_published_articles`, plus category, funder, publisher, and usage statistics) | Date and category browse, one preprint, explicit published DOI | Matches the bioRxiv API, which is an interval browse, and keeps the preprint distinct from the journal article | Nine tools include funder, publisher, and download statistics that are not part of finding or reading a paper. No keyword search. It is a second server the model has to choose |
| OpenAI Codex life-science skills (`ncbi-entrez-skill`, `biorxiv-skill`, `ncbi-pmc-skill` in `openai/plugins`) | Not MCP tools. Compact E-utilities and bioRxiv `details` / `pubs` calls, `max_items` around 10, markdown unless raw output is requested | Small pages and a hard split between preprint metadata and the published DOI | The model has to know API paths. There is no shared publication record. PMC open access is a third skill |
| `@cyanheads/pubmed-mcp-server` (Apache-2.0; design reference only) | Eleven tools: PubMed search and fetch, full text with a PMC then Europe PMC then Unpaywall cascade, Europe PMC search, related `similar` / `cited_by` / `references`, batch citation lookup, id conversion, plus spell check, MeSH, and citation styles | The full-text cascade, truncated search abstracts, named provider, and batch citation lookup are the right agent contracts | The extra tools (spell check, MeSH browser, APA/MLA) are a second product. Europe PMC is a parallel search tool instead of a route |
| Europe PMC MCP servers (for example CoChatAI and pyEuropePMC) | Search, one article, citations, references, JATS full text | Cursor search and a warning that full-text licenses differ by article | One source. Text-mined annotations are a different task from reading a paper |
| OpenAlex and Semantic Scholar MCP servers | `search_works` / `get_work` / citations, often plus authors, institutions, and recommendations | DOI-only works and a real citation graph. Title match is useful | Raw filter grammar and author or institution tools are a poor default for a biomedical reading workflow. Semantic Scholar's unauthenticated rate limit is easy to trip |

`phi-literature` keeps the PubMed search syntax, the citation lookup, the separate preprint and published records, the clipped pages, and the license gate. It drops id conversion, copyright, spell check, MeSH browse, citation formatting, funder search, and usage statistics as their own tools, because those fields are either on the publication or outside the reading tasks above.

## What was not implemented

- **DepMap on the cancer server.** The public models API returns an HTML bot check, so cancer omics does not call it. DepMap release files are handled by the pharmacogenomics pipeline via Figshare, not by a live portal client.
- **ChIP-Atlas.** `https://chip-atlas.dbcls.jp/api/v1/experiment` returned AccessDenied, and the experiment-list keys that were tried were missing. Roadmap Epigenomics is queried through ENCODE instead.
- **CELLxGENE Census matrices and controlled GDC/ENCODE/HCA/HuBMAP/MetaboLights files.** Responses stay metadata. Controlled files have `url: null`.
- **A MetaboLights metabolite catalog.** `/studies/{id}/metabolites` is not a JSON table. The `metabolite` schema exists for a later source.
- **Catalog manifests.** Added under `resources/connectors/` for the five omics servers after the HTTPS gateway went live on `phi.omicsagent.cn`. Pharmacogenomics and literature stay out of the catalog until `/pharmacogenomics/mcp` and `/literature/mcp` are on that host. The connector count stays at 27.
- **Literature PDF parsing, paywalled text, and paid APIs.** `get_fulltext` does not download or parse PDFs. Closed and unknown licenses return links. Semantic Scholar and arXiv are opt-in. OpenCitations is only the fallback for citations and references. bioRxiv and medRxiv keyword search is not invented on top of their date API. A category filter on a bioRxiv date browse applies to the fetched page.

## Planned

HMDB will be added to the metabolomics server from a local bulk download. The public site is behind Cloudflare, so per-metabolite HTTP lookups are not part of this server.

## Add them to the Phi catalog

After the gateway is reachable over HTTPS, add one connector directory per server. Do not point a manifest at a URL that is not public yet. For the deployed host `https://phi.omicsagent.cn`, the transcriptomics manifest is:

```yaml
schemaVersion: 1
id: transcriptomics-single-cell
type: mcp
version: 1.0.0
title: Transcriptomics and single-cell
summary: 检索 GTEx、CELLxGENE、HCA 与 HuBMAP 的开放组学元数据
connector:
  transport: http
  publisher: Phi Packages
  category: 健康与生命科学
  homepage: https://github.com/dxsbiocc/phi-packages
  url: https://phi.omicsagent.cn/transcriptomics/mcp
  auth: none
```

Use the same shape for `cancer-omics` (`/cancer/mcp`), `epigenomics` (`/epigenomics/mcp`), `proteomics` (`/proteomics/mcp`), and `metabolomics` (`/metabolomics/mcp`). Add `pharmacogenomics` (`/pharmacogenomics/mcp`) and `literature` (`/literature/mcp`) only after those paths respond on the public host. Add an original icon that passes Phi's SVG checks, record it outside the 18-icon migration in `docs/content-icons.md` and `SOURCE.json` `addedContent`, and set `SOURCE.json` `counts.connectors` plus both README connector counts to the new directory total. Leave `initialImportCounts.connectors` at 18. Re-run `node --test tests/content-icons.test.mjs tests/content-boundaries.test.mjs`.
