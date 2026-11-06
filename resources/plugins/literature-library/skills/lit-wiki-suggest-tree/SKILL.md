---
name: lit-wiki-suggest-tree
description: >-
  Propose a theme tree from publication tags/titles. Outputs JSON proposals only;
  does not write the vault. Apply with lit-wiki-apply-tree --confirm.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: wiki_suggest_tree
      description: Emit theme-branch proposals JSON (no vault write).
      run: [python, ./scripts/wiki_suggest_tree.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          max_branches:
            type: integer
          out_path:
            type: string
            format: project-path
        required:
          - project_root
          - collection_id
      approval: read
      output: ./schemas/wiki-suggest-tree-result.json
      timeoutSeconds: 60
---

# lit-wiki-suggest-tree

工具名：`lit_wiki_suggest_tree`。**只出提案，不写库。**
