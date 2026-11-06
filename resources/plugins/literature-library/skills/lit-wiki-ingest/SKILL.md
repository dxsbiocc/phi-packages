---
name: lit-wiki-ingest
description: >-
  Ingest a publication into the wiki vault. Wiki prose is produced by the
  invoking agent's model (request → JSON → apply). Cross-paper concept deep
  merge via concept_merge_request / concept_merge_apply. Scripts validate + write.
  No plugin LLM API keys.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
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
---

# lit-wiki-ingest (v0.5.3 — agent compile + concept deep merge + status)

**Wiki prose comes from the invoking agent's already-running model.** This skill never requires a plugin API key.

## Agent steps

```bash
# 1) Emit grounded request (prompt + schema + title/abstract)
python scripts/wiki_compile_request.py \
  --project_root ./demo-project --collection_id col-… \
  --publication_id 'doi:10.1016/j.cell.2017.05.035' --ensure_vault
# → compile_request.json + exit waiting_for_agent

# 2) YOU (the agent) synthesize compile_response.json with your model
#    using only the request's title/abstract/metadata.

# 3) Apply
python scripts/wiki_compile_apply.py \
  --project_root ./demo-project --collection_id col-… \
  --publication_id 'doi:…' --response path/to/compile_response.json
```

Equivalent one-tool path: `wiki_ingest.py --compile agent` (request) then
`wiki_ingest.py --compile agent --llm-json compile_response.json` (apply).

Offline smoke: `--compile mock`. Emergency: `--compile heuristic`.


## Concept deep merge (v0.4.5)

Each wiki ingest writes `wiki/concepts/_contrib/<slug>/<safe_id>.json`. When a
concept already had another source, ingest returns `concept_merge_pending` and a
soft stub under `raw/concept_merge_request_<slug>.json` (no invented synthesis).

```bash
python scripts/concept_merge_request.py --project_root P --collection_id C --slug t-cell-exhaustion
# agent fills concept_merge_response_*.json
python scripts/concept_merge_apply.py --project_root P --collection_id C \
  --slug t-cell-exhaustion --llm-json path/to/response.json
# offline:
python scripts/concept_merge_apply.py ... --slug t-cell-exhaustion --mock
```

Rewrites Definition + ## Across papers; preserves ## Sources. `needs_review`
contributions require `needs_review_banner`. Optional `canonical_slug` / `merge_into`
→ redirect note on old page.

## Agent status helper (0.5.3)

```bash
python scripts/wiki_agent_status.py --project_root P --collection_id C
python scripts/wiki_agent_status.py --project_root P --collection_id C \
  --emit_requests --limit 5 --out_dir ./wiki-requests
```

Lists enriched pubs missing `wiki/sources/*` and concepts with merge pending / ≥2 contribs.
`--emit_requests` only writes compile_request.json — **you** still fill responses (no API keys).

`related_concepts` on merge apply resolve to existing concept slugs (aliases like `nsclc`→`non-small-cell-lung-cancer`); unresolved dropped + reported as `related_unresolved`.

工具名：`lit_wiki_ingest` / `lit_wiki_compile_request` / `lit_wiki_compile_apply` /
`lit_wiki_agent_status` / `lit_concept_merge_request` / `lit_concept_merge_apply`。
