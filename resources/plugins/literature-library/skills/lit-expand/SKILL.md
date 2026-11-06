---
name: lit-expand
description: >-
  Seed → related / citations / references via in-plugin HTTP (PubMed elink,
  OpenAlex related/referenced/cites, OpenCitations fallback). Preview by default;
  ingest only with --confirm_write. No MCP server; no plugin API keys.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: expand
      description: >-
        Expand from seed_id. Default prints hits JSON without DB write;
        --confirm_write writes expands/<id>/hits.jsonl and ingests.
      run: [python, ./scripts/expand.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          seed_id:
            type: string
          relation:
            type: string
            enum: [similar, citations, references]
          limit:
            type: integer
          confirm_write:
            type: boolean
          dry_run:
            type: boolean
          mock:
            type: boolean
          skip_ingest:
            type: boolean
          no_enrich:
            type: boolean
          max_rounds:
            type: integer
          low_yield_ratio:
            type: number
          seeds_per_round:
            type: integer
        required:
          - project_root
          - collection_id
          - seed_id
      approval: read
      output: ./schemas/expand-result.json
      timeoutSeconds: 120
---

# lit-expand

从一篇种子文献扩展 **similar / citations / references**。草案走插件内 HTTP（`providers_http`），**不**内嵌 MCP、**不**打包 API key。

## 行为

| 参数 | 说明 |
|--|--|
| `--seed_id` | DOI / PMID / PMCID / arXiv / OpenAlex `W…` / 库内 publication id |
| `--relation` | `similar`（默认）/ `citations` / `references` |
| `--limit` | 1–50，默认 10 |
| 默认 | **只预览** hits + meta；`wrote: false` |
| `--confirm_write` | 写入 `collections/<id>/expands/<exp-id>/hits.jsonl` 并 ingest |
| `--dry-run` | 强制预览（覆盖 confirm） |
| `--mock` | 离线假数据 |
| `--skip_ingest` | 只写 jsonl，不入库 |
| `--no_enrich` | 入库后不跑 lit-enrich |
| `--max_rounds` | 雪球轮次（默认 1）；触顶未饱和 → `max_rounds_not_saturated` |
| `--low_yield_ratio` | 低产阈值（默认 0.05）；连续 2 轮低产 → saturated |
| `--seeds_per_round` | 每轮最多扩展的种子数（默认 5） |

## 关系 → 后端

| relation | 主路径 | 回退 |
|--|--|--|
| similar | PMID → PubMed elink `pubmed_pubmed` | OpenAlex `related_works` hydrate |
| references | OpenAlex `referenced_works` hydrate | OpenCitations references → OpenAlex DOI hydrate |
| citations | OpenAlex filter `cites:W…` | OpenCitations citations → OpenAlex DOI hydrate |

工具名：`lit_expand`。stdout `stub: false`（happy path）。

输出含 `saturated` / `stop_reason` / `rounds`。**禁止**在未满足饱和判据时宣称闭包。
