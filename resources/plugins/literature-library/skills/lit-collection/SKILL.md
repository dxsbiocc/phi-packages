---
name: lit-collection
description: >-
  Create or list literature collections (课题) inside a project. Initializes
  library.sqlite, an empty theme tree root, and (by default) an LLM Wiki vault
  under collections/<id>/wiki-vault/. Pass --no_wiki to skip vault init.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: collection
      description: Create or list collections under a project root (library.sqlite).
      run: [python, ./scripts/collection.py]
      args:
        type: object
        additionalProperties: false
        properties:
          action:
            type: string
            enum: [create, list]
          project_root:
            type: string
            format: project-path
            description: Project directory that will hold library.sqlite
          name:
            type: string
            description: Collection name (required for create)
          description:
            type: string
          no_wiki:
            type: boolean
            description: If true, skip wiki-vault init on create
        required:
          - action
          - project_root
      approval: write
      output: ./schemas/collection-result.json
      timeoutSeconds: 60
---

# lit-collection

课题 = collection。每个 Phi 项目一份 `library.sqlite`；collection 目录下挂 `runs/`、`files/`、默认 `wiki-vault/`。

```bash
python scripts/collection.py --action create --project_root ./demo-project --name "HCC T cells"
python scripts/collection.py --action create --project_root ./demo-project --name "No wiki" --no_wiki
python scripts/collection.py --action list --project_root ./demo-project
```

工具名：`lit_collection`。
