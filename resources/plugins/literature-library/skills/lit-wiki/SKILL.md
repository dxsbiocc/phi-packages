---
name: lit-wiki
description: >-
  Compile one collection's LLM Wiki (raw snapshots, wiki pages, AGENTS.md, theme tree). The invoking model writes summaries, claims, and cross-paper concept merges; these tools build the grounded prompt, validate JSON, and write the vault. Use for wiki init, compile request/apply, concept deep-merge, compile status, keyword query, lint, and theme-tree suggest then confirm. Paper search and SQLite ingest belong to the literature-library skill. No embeddings and no plugin API keys.
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
    - name: wiki_agent_status
      description: >-
        List enriched pubs missing wiki sources and concepts awaiting merge;
        optional batch emit compile_request.json (agent fills responses).
      run: [python, ./scripts/wiki_agent_status.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          limit:
            type: integer
          emit_requests:
            type: boolean
          out_dir:
            type: string
            format: project-path
        required:
          - project_root
          - collection_id
      approval: read
      timeoutSeconds: 120
    - name: wiki_compile_request
      description: Emit grounded compile_request.json for the invoking agent to fill.
      run: [python, ./scripts/wiki_compile_request.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          publication_id:
            type: string
          ensure_vault:
            type: boolean
          emit_request:
            type: string
            format: project-path
            description: "Also accepted as --request-out"
        required:
          - project_root
          - collection_id
          - publication_id
      approval: write
      timeoutSeconds: 60
    - name: wiki_compile_apply
      description: Validate agent compile JSON and write wiki source/concept/entity pages.
      run: [python, ./scripts/wiki_compile_apply.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          publication_id:
            type: string
          response:
            type: string
            format: project-path
          llm_json:
            type: string
            format: project-path
          apply_compile:
            type: string
            format: project-path
          allow_heuristic_fallback:
            type: boolean
        required:
          - project_root
          - collection_id
          - publication_id
      approval: write
      timeoutSeconds: 120
    - name: wiki_ingest
      description: >-
        Ingest one publication. Default --compile agent: emit compile_request
        if no --llm-json; apply when agent JSON is provided.
      run: [python, ./scripts/wiki_ingest.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          publication_id:
            type: string
          raw_path:
            type: string
            format: project-path
          ensure_vault:
            type: boolean
          compile:
            type: string
            description: "agent (default) | mock | heuristic | file"
          llm_json:
            type: string
            format: project-path
            description: Agent-produced compile_response.json
          apply_compile:
            type: string
            format: project-path
          request_only:
            type: boolean
          allow_heuristic_fallback:
            type: boolean
        required:
          - project_root
          - collection_id
      approval: write
      output: ./schemas/wiki-ingest-result.json
      timeoutSeconds: 120
    - name: concept_merge_request
      description: >-
        Emit grounded concept-merge request for concepts with ≥2 per-source
        contributions (agent fills JSON; no plugin API keys).
      run: [python, ./scripts/concept_merge_request.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          vault:
            type: string
            format: project-path
          slug:
            type: string
          all:
            type: boolean
          emit_request:
            type: string
            format: project-path
        required: []
      approval: write
      output: ./schemas/concept-merge-result.json
      timeoutSeconds: 60
    - name: concept_merge_apply
      description: Validate agent/mock concept-merge JSON and rewrite concept page.
      run: [python, ./scripts/concept_merge_apply.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          vault:
            type: string
            format: project-path
          slug:
            type: string
          llm_json:
            type: string
            format: project-path
          mock:
            type: boolean
        required:
          - slug
      approval: write
      output: ./schemas/concept-merge-result.json
      timeoutSeconds: 60
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

# LLM Wiki

Compile layer for one collection vault at `collections/<id>/wiki-vault/` (`raw/`, `wiki/`, `AGENTS.md`, `wiki/index.md`, `wiki/log.md`). Raw snapshots stay immutable. Wiki pages are what a reader should open.

Search, file ingest, and shell enrich stay on the `literature-library` skill (`lit_search`, `lit_ingest`, `lit_enrich`). This skill does not embed vectors and does not hold an API key.

The model that is already running writes Summary, key claims, and cross-paper Definition / Across papers. These tools build a grounded prompt and schema, validate JSON, and write the vault. Do not set `OPENAI_API_KEY` or `LIT_WIKI_LLM_API_KEY` for this plugin.

The old `lit-wiki-*` skill folders are this skill. Tool ids are unchanged: `lit_wiki_init`, `lit_wiki_ingest`, `lit_wiki_compile_request`, `lit_wiki_compile_apply`, `lit_concept_merge_request`, `lit_concept_merge_apply`, `lit_wiki_agent_status`, `lit_wiki_query`, `lit_wiki_lint`, `lit_wiki_suggest_tree`, `lit_wiki_apply_tree`.

## Compile one paper

If the row is only a DOI or PMID, the title is a placeholder, or the abstract is empty, run `lit_enrich` first. A request that carries `warnings` or `prefer_enrich_first` means the same thing.

1. `lit_wiki_compile_request`, or `lit_wiki_ingest` without `llm_json`, writes `compile_request.json` and returns `waiting_for_agent`. `ensure_vault` creates the skeleton when it is missing.
2. Read that request's title, abstract, and metadata. Write `compile_response.json` with the current model, using only those fields. Do not invent PMIDs or statistics. Do not emit fragment symbols `l1` or `ifn` from PD-L1 or IFN-γ. Map layilin to `LAYN`.
3. `lit_wiki_compile_apply`, or `lit_wiki_ingest` with `llm_json`, writes `wiki/sources`, concepts, and entities.

`lit_wiki_init` creates the same skeleton (idempotent). `lit_collection` create on the library skill already calls it unless `no_wiki` was set.

Offline smoke only: `compile=mock`. Emergency scaffolding, not target quality: `compile=heuristic` or `LIT_WIKI_COMPILE=heuristic`. Do not recommend mock or heuristic for a live vault.

`lit_wiki_agent_status` lists enriched publications that still have no `wiki/sources` page, and concepts that have a merge pending or at least two contributions. `emit_requests` only writes `compile_request.json`. You still fill the responses.

## Concept deep merge

One paper's compile upserts links and writes `wiki/concepts/_contrib/<slug>/<safe_id>.json`. It does not synthesize Definition or Across papers. When ingest returns `concept_merge_pending` (and a stub under `raw/concept_merge_request_<slug>.json`):

1. `lit_concept_merge_request` for that slug, or `all`.
2. Fill the merge JSON from the contributions plus the current Definition and Sources. Do not invent synthesis in a script.
3. `lit_concept_merge_apply` with `llm_json`. Offline: `mock`.

The apply rewrites Definition and `## Across papers` and keeps `## Sources`. A `needs_review` source must keep `needs_review_banner`. Optional `canonical_slug` / `merge_into` leaves a redirect on the old page. `related_concepts` must resolve to an existing slug (an alias such as `nsclc` can resolve); unresolved names are dropped and reported as `related_unresolved`.

## Query, lint, and the theme tree

`lit_wiki_query` is keyword AND across wiki markdown plus SQLite FTS or LIKE. `wiki_only` skips the library. There is no vector retrieval.

`lit_wiki_lint` reports orphan pages, index gaps, broken relative links, stale `needs_review` banners, and noisy entities.

`lit_wiki_suggest_tree` emits branch proposals from tags and titles and does not write the vault. Read the proposal with the user. `lit_wiki_apply_tree` is a dry-run (`wrote: false`) until `confirm`. Only the confirmed call updates `wiki/index.md` and `theme_nodes`.

## Repoint

On a publication-id upgrade or a version merge, `wiki_repoint_publication` renames `wiki/sources/<old_safe>.md` (and the raw file and contrib), updates index and log links, and leaves a stub at the old path. `lit_versions` apply/confirm and `lit_merge_ids` apply call this unless repoint is turned off. Do not hand-rename those pages.

## Shared modules

`scripts/wiki_lib.py` and `scripts/wiki_llm.py` are the vault layout, lint, tree, and JSON-validation helpers. The tools in this skill import them. They are not separate tools.
