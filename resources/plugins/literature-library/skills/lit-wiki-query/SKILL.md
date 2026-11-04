---
name: lit-wiki-query
description: >-
  Keyword search over wiki markdown plus SQLite FTS/LIKE metadata. No embeddings
  or vector retrieval.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: wiki_query
      description: Search wiki pages and optional SQLite publications by keywords.
      run: [python, ./scripts/wiki_query.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          q:
            type: string
            description: "Keyword query (AND terms); CLI also accepts --query"
          query:
            type: string
            description: "Alias of q"
          limit:
            type: integer
          wiki_only:
            type: boolean
        required:
          - project_root
          - collection_id
          - q
      approval: read
      output: ./schemas/wiki-query-result.json
      timeoutSeconds: 60
---

# lit-wiki-query

工具名：`lit_wiki_query`。**无向量检索。** 多词为 **AND** 语义。

```bash
python scripts/wiki_query.py \
  --project_root ./demo-project --collection_id col-... \
  --query "exhaustion LAYN"
```
