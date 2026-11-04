---
name: lit-wiki-lib
description: >-
  Shared LLM-Wiki helpers (vault layout, ingest templates, lint, suggest/apply
  theme tree). Other lit-wiki-* skills import scripts/wiki_lib.py. No embeddings.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts: []
---

# lit-wiki-lib

内部共享：`scripts/wiki_lib.py`。Vault = `collections/<id>/wiki-vault/{raw,wiki,AGENTS.md}`。

**无 embedding。** TextIn 文档树不是知识库骨架。

## v0.2.2 entities

Fragment guard: no `l1`/`ifn` from PD-L1/IFN-γ; curated allowlist + `noisy_entity` lint.

## v0.4.5 concept merge

`wiki_lib.py`: per-source contributions under `wiki/concepts/_contrib/`.
`wiki_llm.py`: `build_concept_merge_request` / `mock_concept_merge_payload` / `apply_concept_merge`.

## v0.4.8 wiki_repoint

`wiki_repoint_publication` / `wiki_repoint_remaps`: on publication id upgrade or version merge, rename `wiki/sources/<old_safe>.md` (+ raw + contrib), update index/log links, leave stub at old path. Called from lit-versions apply/confirm and lit-db merge_ids --apply.
