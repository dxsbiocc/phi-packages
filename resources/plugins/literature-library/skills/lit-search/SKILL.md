---
name: lit-search
description: >-
  Execute a configured literature run: fetch hits via in-plugin HTTP adapters
  (pubmed, europepmc, openalex, crossref, arxiv, biorxiv, medrxiv, pmc, doaj,
  datacite, zenodo; semantic_scholar soft-fails without S2_API_KEY) or --mock;
  write hits.jsonl; auto-ingest with dedup. Draft does not ship an MCP server;
  production Phi may later swap to MCP search_publications (phi-literature).
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: search
      description: >-
        Run search for an existing run_id; write hits.jsonl; auto-call ingest.
        Use --mock for offline smoke. Optional --sources overrides run.json.
      run: [python, ./scripts/search.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          run_id:
            type: string
          mock:
            type: boolean
            description: Offline mock hits (no network)
          skip_ingest:
            type: boolean
          sources:
            type: string
            description: Comma-separated source override (e.g. pmc,doaj,datacite,zenodo)
        required:
          - project_root
          - collection_id
          - run_id
      approval: write
      output: ./schemas/search-result.json
      timeoutSeconds: 180
---

# lit-search

1. 读 `runs/<run_id>/run.json`（可用 `--sources` 覆盖）
2. 拉取 hits → `hits.jsonl`（`providers_http.py` 公共 REST）
3. **自动**调用 ingest 去重入库（`--skip_ingest` 可关）

## Provider 状态（0.4.9-draft）

| 源 | 本草案（插件内 HTTP） | 说明 |
|--|--|--|
| pubmed | **live** E-utilities | |
| europepmc | **live** REST | |
| openalex | **live** api.openalex.org | 可选 `OPENALEX_MAILTO` / `CONTACT_EMAIL` |
| crossref | **live** api.crossref.org | 同上 mailto |
| arxiv | **live** export.arxiv.org Atom | 默认偏 q-bio；无命中再 `all:` |
| biorxiv / medrxiv | **live** | 无关键词 API：DOI detail / 日期窗 / Europe PMC `JOURNAL` 回退 |
| pmc | **live** E-utilities `db=pmc` | esearch + esummary；映射 PMCID/DOI/PMID |
| doaj | **live** doaj.org/api/search/articles | OA 期刊文章；形状变化 soft-fail |
| datacite | **live** api.datacite.org/dois | 数据集/预印本等 DOI；礼貌 mailto |
| zenodo | **live** zenodo.org/api/records | 论文/预印本/数据集 |
| semantic_scholar / S2 | **live or soft-fail** | 需环境变量 `S2_API_KEY`（不入库）；401/403/429 → 空结果 + meta |
| `--mock` | 离线 2 条假数据 | 测试用 |

不接 CORE / Scopus / WoS / SmartLib / CNKI 等商业网关。

**架构：** 对外只暴露 `lit-search` CLI；**不**在插件内嵌第二套 MCP server。生产 Phi *可以* 以后把同一 CLI 换成 MCP client 调 `phi-literature`，但草案路径是 HTTP。

工具名：`lit_search`。

## 0.5.0

Optional `--close_run` / `--min_unique_pubs` / `--fail_under_coverage` after ingest; writes durable `hits.meta.json`.


## Backends (0.6.0)

| Backend | How |
|--|--|
| `http` (default) | `providers_http.py` |
| `mcp` | Optional client `providers_mcp.py` → phi-literature `search_publications` via `LIT_MCP_CMD` or `--mcp_mock_fixture` |
| env | `LIT_SEARCH_BACKEND=mcp` |

This plugin does **not** ship an MCP server. Soft-fail / `--mcp_fallback_http` when MCP is unavailable.
