---
name: lit-query
description: >-
  Query project library.sqlite by collection, FTS/LIKE on title+abstract+tags,
  optional tag and needs_review filters. No embedding in v1.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: query
      description: Search publications in library.sqlite; JSON stdout.
      run: [python, ./scripts/query.py]
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
            description: "FTS/LIKE query; CLI also accepts --query"
          query:
            type: string
            description: "Alias of q"
          tag:
            type: string
          needs_review:
            type: string
            description: "0 or 1"
          limit:
            type: integer
        required:
          - project_root
      approval: read
      output: ./schemas/query-result.json
      timeoutSeconds: 60
---

# lit-query

v1：**SQL + FTS5**（不可用则 `LIKE`）。无 embedding。

```bash
python scripts/query.py --project_root ./demo-project --collection_id col-... --query "exhaustion"
```

工具名：`lit_query`。
