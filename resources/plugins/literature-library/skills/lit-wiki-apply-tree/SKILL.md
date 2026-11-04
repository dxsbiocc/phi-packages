---
name: lit-wiki-apply-tree
description: >-
  Apply theme-tree proposals to wiki/index.md and theme_nodes only when
  --confirm is set. Without --confirm, dry-run (wrote=false).
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: wiki_apply_tree
      description: Apply proposals JSON to wiki index; requires --confirm to write.
      run: [python, ./scripts/wiki_apply_tree.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          proposals_path:
            type: string
            format: project-path
          confirm:
            type: boolean
        required:
          - project_root
          - collection_id
          - proposals_path
      approval: write
      output: ./schemas/wiki-apply-tree-result.json
      timeoutSeconds: 120
---

# lit-wiki-apply-tree

```bash
# dry-run
python scripts/wiki_apply_tree.py --project_root ./demo --collection_id col-… \
  --proposals_path ./proposals.json
# write
python scripts/wiki_apply_tree.py … --proposals_path ./proposals.json --confirm
```

工具名：`lit_wiki_apply_tree`。
