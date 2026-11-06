---
name: lit-wiki-init
description: >-
  Create LLM Wiki vault under a collection: wiki-vault/{raw,wiki,AGENTS.md,
  wiki/index.md, wiki/log.md}. Collection ≈ one vault root. No embeddings.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: wiki_init
      description: Initialize wiki-vault skeleton for a collection (idempotent).
      run: [python, ./scripts/wiki_init.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
        required:
          - project_root
          - collection_id
      approval: write
      output: ./schemas/wiki-init-result.json
      timeoutSeconds: 60
---

# lit-wiki-init

```bash
python scripts/wiki_init.py --project_root ./demo-project --collection_id col-…
```

工具名：`lit_wiki_init`。
