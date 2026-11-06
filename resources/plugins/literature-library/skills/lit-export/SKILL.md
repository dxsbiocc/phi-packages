---
name: lit-export
description: >-
  Export a collection to search_log.jsonl for research-evidence consumers.
  Does not write hypothesis cards or run scientific-gate.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: export
      description: Write search_log.jsonl for a collection.
      run: [python, ./scripts/export.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          out_path:
            type: string
            format: project-path
        required:
          - project_root
          - collection_id
          - out_path
      approval: write
      output: ./schemas/export-result.json
      timeoutSeconds: 60
---

# lit-export

给 `research-evidence` 消费的 `search_log.jsonl`。本插件**不**写假设卡 / gate。

工具名：`lit_export`。
