# Omics MCP servers

Four read-only MCP servers for omics databases that do not ship an official MCP endpoint. They share one entity model (`docs/schemas.md`) and one HTTP client (timeouts, retries, per-host rate limit, on-disk cache).

| Server | Package | Sources |
| --- | --- | --- |
| Transcriptomics and single-cell | `phi-omics-transcriptomics` | GTEx v10, CELLxGENE Discover, CellGuide, HCA Azul, HuBMAP search |
| Cancer omics | `phi-omics-cancer` | GDC, cBioPortal public REST, PDC GraphQL |
| Epigenomics | `phi-omics-epigenomics` | ENCODE portal, Roadmap via `award.rfa=Roadmap` |
| Proteomics and metabolomics | `phi-omics-proteomics` | PRIDE, Human Protein Atlas, MetaboLights |

The same TypeScript runs locally over stdio or Streamable HTTP, and in Docker on one host. Packages are private workspace members. An npm scope and a publish workflow are intentionally not set up yet; `bin`, `files`, and `engines` are already in each package so that step can be added later without moving code.

These servers are not catalog connectors. `resources/connectors/` is unchanged until a public HTTPS URL exists. See [Add them to the Phi catalog](#add-them-to-the-phi-catalog) after deployment.

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
```

The bins are `phi-omics-transcriptomics`, `phi-omics-cancer`, `phi-omics-epigenomics`, `phi-omics-proteomics`, and `phi-omics-gateway`. Bun can run the same `dist/cli.js` files; there are no native addons.

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
| `PHI_OMICS_MAX_BODY_BYTES` | `8000000` | Reject larger upstream bodies |
| `GTEX_DATASET` | `gtex_v10` | GTEx dataset |
| `GTEX_GENCODE` | `v39` | GTEx GENCODE version |
| `GTEX_GENOME` | `GRCh38/hg38` | GTEx genome build argument |
| `CELLGUIDE_SNAPSHOT` | `1692975365` | Pinned CellGuide public snapshot |
| `MCP_SERVER` | `all` | Gateway only: one server name, or `all` |

No upstream API key is required. The servers do not send tokens and do not download controlled files.

## One host with Docker

From `mcp-servers/`:

```bash
docker compose up --build
```

The image listens on port 8080:

- `GET /health`
- `POST /transcriptomics/mcp`
- `POST /cancer/mcp`
- `POST /epigenomics/mcp`
- `POST /proteomics/mcp`

`MCP_SERVER=transcriptomics-single-cell` (or `cancer-omics`, `epigenomics`, `proteomics-metabolomics`) serves that one server at `/mcp` instead.

Put TLS on a reverse proxy. `Caddyfile.example` proxies a hostname to `127.0.0.1:8080`. Caddy obtains and renews certificates. Nginx or another proxy is the same pattern: forward `POST`, `GET`, and `DELETE` to the Node process, preserve `Authorization` and `Mcp-Session-Id` if you add auth later, and do not buffer streaming responses. This process does not terminate TLS and does not authenticate clients. Bind it to localhost when the proxy is on the same machine.

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

### Proteomics and metabolomics

`usage_terms`, `pride_search_projects`, `pride_get_project`, `pride_list_files`, `hpa_search_genes`, `hpa_get_gene`, `hpa_list_proteins`, `metabolights_list_studies`, `metabolights_get_study`, `metabolights_list_files`, `metabolights_list_assays`.

## What was not implemented

- **DepMap.** `https://depmap.org/portal/api/public/models` returned an HTML bot-check page, not JSON.
- **ChIP-Atlas.** `https://chip-atlas.dbcls.jp/api/v1/experiment` returned AccessDenied, and the experiment-list keys that were tried were missing. Roadmap Epigenomics is queried through ENCODE instead.
- **CELLxGENE Census matrices and controlled GDC/ENCODE/HCA/HuBMAP/MetaboLights files.** Responses stay metadata. Controlled files have `url: null`.
- **A MetaboLights metabolite catalog.** `/studies/{id}/metabolites` is not a JSON table. The `metabolite` schema exists for a later source.
- **Catalog manifests.** There is no public URL yet, so nothing was added under `resources/connectors/`.

## Add them to the Phi catalog

After the gateway is reachable over HTTPS, add one connector directory per server. Do not point a manifest at a URL that is not public yet. For a host `https://omics.example.com`, the transcriptomics manifest is:

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
  url: https://omics.example.com/transcriptomics/mcp
  auth: none
```

Use the same shape for `cancer-omics` (`/cancer/mcp`), `epigenomics` (`/epigenomics/mcp`), and `proteomics-metabolomics` (`/proteomics/mcp`). Add an original icon that passes Phi's SVG checks, record it outside the 18-icon migration in `docs/content-icons.md` and `SOURCE.json` `addedContent`, and set `SOURCE.json` `counts.connectors` plus both README connector counts to the new directory total. Leave `initialImportCounts.connectors` at 18. Re-run `node --test tests/content-icons.test.mjs tests/content-boundaries.test.mjs`.
