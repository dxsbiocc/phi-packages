---
name: lit-run-config
description: >-
  Freeze a literature search run protocol into runs/<run_id>/run.yaml + run.json
  (query, sources, dates, limit_per_source, provider_notes, frozen_at). Open HTTP
  sources only (draft ≤8; includes pmc/doaj/datacite/zenodo). No commercial sources.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: run_config
      description: Freeze run protocol (query, sources, dates, type, limit).
      run: [python, ./scripts/run_config.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          query:
            type: string
          sources:
            type: string
            description: Comma-separated sources (max 8; e.g. pubmed,pmc,doaj,datacite)
          date_from:
            type: string
          date_to:
            type: string
          type:
            type: string
          limit:
            type: integer
          run_id:
            type: string
          gold_set:
            type: string
            description: Comma-separated strong IDs or @file; guidance 5–10
          out_of_window_controls:
            type: string
            description: Comma-separated 1–2 out-of-window control IDs or @file
        required:
          - project_root
          - collection_id
          - query
      approval: write
      output: ./schemas/run-config-result.json
      timeoutSeconds: 60
---

# lit-run-config

一次检索 = 一个 run。写入 `collections/<id>/runs/<run_id>/run.yaml` + **冻结** `run.json`：

- `protocol_version` (0.5.2), `query`, `sources`, `date_from`/`date_to`
- `limit` / `limit_per_source`, `provider_notes`, `frozen_at`
- `gold_set`, `out_of_window_controls` (caller verifies gold DOIs; plugin format-checks only)
- `no_commercial_sources: true`, `status: configured`

允许源：pubmed, europepmc, openalex, crossref, semantic_scholar, arxiv, biorxiv, medrxiv, pmc, doaj, datacite, zenodo（最多 8）。

工具名：`lit_run_config`。
