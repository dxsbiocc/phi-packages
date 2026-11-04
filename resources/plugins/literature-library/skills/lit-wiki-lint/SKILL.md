---
name: lit-wiki-lint
description: >-
  Lint wiki vault for orphan pages, missing index entries, broken relative
  links, and stale needs_review banners.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: wiki_lint
      description: Lint a collection wiki-vault and report issues JSON.
      run: [python, ./scripts/wiki_lint.py]
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
      approval: read
      output: ./schemas/wiki-lint-result.json
      timeoutSeconds: 120
---

# lit-wiki-lint

工具名：`lit_wiki_lint`。
