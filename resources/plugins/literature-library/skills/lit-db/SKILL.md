---
name: lit-db
description: >-
  Shared SQLite helpers for literature-library (init_db, normalize_*,
  resolve_or_create, query_publications). Not usually called alone —
  lit-collection / lit-ingest / lit-query import scripts/db.py.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: init_db
      description: Initialize library.sqlite schema (idempotent). Returns path and FTS flag.
      run: [python, ./scripts/init_db_cli.py]
      args:
        type: object
        additionalProperties: false
        properties:
          db_path:
            type: string
            format: project-path
            description: Path to library.sqlite
        required:
          - db_path
      approval: write
      output: ./schemas/init-db-result.json
      timeoutSeconds: 30
    - name: merge_ids
      description: >-
        Dry-run (default) or --apply high-confidence merges of complementary
        DOI/PMID publication clusters into DOI-preferred canonical ids.
      run: [python, ./scripts/merge_ids_cli.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          db_path:
            type: string
            format: project-path
          apply:
            type: boolean
            description: Write merges when true; default dry-run
          include_low_confidence:
            type: boolean
        required: []
      approval: write
      timeoutSeconds: 60
---

# lit-db

内部共享库：`scripts/db.py`。其他技能通过 `sys.path` 引入。

## 表

collections · publications（含 oa_status/oa_json）· pub_files · aliases · collection_pubs · runs · run_hits · theme_nodes · theme_pubs · publication_versions · pubs_fts (FTS5，不可用则 LIKE)

## 去重

强 ID：DOI > PMID > PMCID > arXiv > OpenAlex → 命中则 merge。  
软指纹（title+first_author+year）→ `needs_review`，**不**自动合并。  
真冲突强 ID → `needs_review`。互补强 ID 空壳/同篇 → 自动合并（见 `merge_ids`）。

