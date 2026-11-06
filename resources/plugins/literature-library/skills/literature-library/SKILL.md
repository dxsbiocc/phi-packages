---
name: literature-library
description: >-
  Find, import, and manage one project literature library. Create a collection, freeze a search run, search (in-plugin HTTP or mock; optional MCP client), ingest files, deduplicate, enrich ID shells, expand from a seed, query, export search_log, report coverage, close a run, list the identity review queue, and link preprint to journal versions. Includes own-account institutional full text and content-hash path rebind. Use for 文献库、课题检索、入库、去重、coverage、版本链. Wiki prose is the lit-wiki skill. No commercial gateways, Google Scholar scrape, embeddings, or plugin API keys.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: init_db
      description: Initialize library.sqlite schema (idempotent). Returns path and FTS flag.
      run: [python, ./scripts/init_db_cli.py]
      args:
        type: object
        additionalProperties: false
        properties:
          db_path:
            type: string
            format: project-path
            description: Path to library.sqlite
        required:
          - db_path
      approval: write
      output: ./schemas/init-db-result.json
      timeoutSeconds: 30
    - name: merge_ids
      description: >-
        Dry-run (default) or --apply high-confidence merges of complementary
        DOI/PMID publication clusters into DOI-preferred canonical ids.
      run: [python, ./scripts/merge_ids_cli.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          db_path:
            type: string
            format: project-path
          apply:
            type: boolean
            description: Write merges when true; default dry-run
          include_low_confidence:
            type: boolean
        required: []
      approval: write
      timeoutSeconds: 60
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
    - name: run_config
      description: Freeze run protocol (query, sources, dates, type, limit).
      run: [python, ./scripts/run_config.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          query:
            type: string
          sources:
            type: string
            description: Comma-separated sources (max 8; e.g. pubmed,pmc,doaj,datacite)
          date_from:
            type: string
          date_to:
            type: string
          type:
            type: string
          limit:
            type: integer
          run_id:
            type: string
          gold_set:
            type: string
            description: Comma-separated strong IDs or @file; guidance 5–10
          out_of_window_controls:
            type: string
            description: Comma-separated 1–2 out-of-window control IDs or @file
        required:
          - project_root
          - collection_id
          - query
      approval: write
      output: ./schemas/run-config-result.json
      timeoutSeconds: 60
    - name: search
      description: >-
        Run search for an existing run_id; write hits.jsonl; auto-call ingest.
        Use --mock for offline smoke. Optional --sources overrides run.json.
      run: [python, ./scripts/search.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          run_id:
            type: string
          mock:
            type: boolean
            description: Offline mock hits (no network)
          skip_ingest:
            type: boolean
          sources:
            type: string
            description: Comma-separated source override (e.g. pmc,doaj,datacite,zenodo)
        required:
          - project_root
          - collection_id
          - run_id
      approval: write
      output: ./schemas/search-result.json
      timeoutSeconds: 180
    - name: ingest
      description: Parse input file and upsert into library.sqlite with dedup.
      run: [python, ./scripts/ingest.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          input_path:
            type: string
            format: input-path
          format:
            type: string
            description: hits|bib|ris|ids|pdf_text|pdf (optional auto-detect)
          run_id:
            type: string
          no_enrich:
            type: boolean
            description: Skip auto-enrich after local/ids ingest (default false)
          enrich_mock:
            type: boolean
            description: Use offline mock fixture when auto-enriching
          verify_oa:
            type: boolean
            description: Soft-verify candidate OA/PDF URLs; never store PDF bytes
          oa_mock_fixture:
            type: string
            format: input-path
            description: JSON map url→mock HEAD result for offline OA verify
        required:
          - project_root
          - collection_id
          - input_path
      approval: write
      output: ./schemas/ingest-result.json
      timeoutSeconds: 120
    - name: rebind_files
      description: >-
        Rebind pub_files.file_path by content_hash after files move (dry-run
        default; --confirm writes).
      run: [python, ./scripts/rebind_files.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          scan_dir:
            type: string
            format: input-path
          confirm:
            type: boolean
          sample_limit:
            type: integer
        required:
          - project_root
          - scan_dir
      approval: write
      timeoutSeconds: 300
    - name: institutional_access
      description: >-
        Mount an already-authorized local PDF, or fetch via LIT_INSTITUTIONAL_PROXY
        / ambient VPN (mockable). Records access_origin + access_method on pub_files.
      run: [python, ./scripts/institutional_access.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          institutional_access:
            type: boolean
            description: Required gate (must be true)
          method:
            type: string
            description: vpn_proxy|session_reuse|saml_browser|user_provided_file
          file_path:
            type: string
            format: input-path
          url:
            type: string
          dest_dir:
            type: string
            format: project-path
          doi:
            type: string
          pmid:
            type: string
          title:
            type: string
          mock_fixture:
            type: string
            format: input-path
          proxy:
            type: string
          no_ocr:
            type: boolean
          ocr_cmd:
            type: string
        required:
          - project_root
          - collection_id
          - institutional_access
          - method
      approval: write
      timeoutSeconds: 300
    - name: enrich
      description: Fetch metadata for shell publications and merge into SQLite.
      run: [python, ./scripts/enrich.py]
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
          dry_run:
            type: boolean
          mock:
            type: boolean
          mock_fixture:
            type: string
            format: project-path
          force:
            type: boolean
          verify_oa:
            type: boolean
            description: Soft-verify candidate OA URLs after enrich; never store PDF bytes
          oa_mock_fixture:
            type: string
            format: project-path
        required:
          - project_root
          - collection_id
      approval: write
      output: ./schemas/enrich-result.json
      timeoutSeconds: 180
    - name: expand
      description: >-
        Expand from seed_id. Default prints hits JSON without DB write;
        --confirm_write writes expands/<id>/hits.jsonl and ingests.
      run: [python, ./scripts/expand.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          seed_id:
            type: string
          relation:
            type: string
            enum: [similar, citations, references]
          limit:
            type: integer
          confirm_write:
            type: boolean
          dry_run:
            type: boolean
          mock:
            type: boolean
          skip_ingest:
            type: boolean
          no_enrich:
            type: boolean
          max_rounds:
            type: integer
          low_yield_ratio:
            type: number
          seeds_per_round:
            type: integer
        required:
          - project_root
          - collection_id
          - seed_id
      approval: read
      output: ./schemas/expand-result.json
      timeoutSeconds: 120
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
    - name: run_close
      description: Coverage gate + mark run status after search/ingest.
      run: [python, ./scripts/run_close.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          db_path:
            type: string
            format: project-path
          run_id:
            type: string
          min_unique_pubs:
            type: integer
          max_needs_review:
            type: integer
          max_version_pending:
            type: integer
          max_unenriched_shells:
            type: integer
          fail_under_coverage:
            type: boolean
          coverage_json:
            type: string
            format: project-path
          coverage_md:
            type: string
            format: project-path
        required:
          - run_id
      approval: write
      output: ./schemas/run-close-result.json
      timeoutSeconds: 60
    - name: coverage
      description: JSON coverage report; optional markdown summary.
      run: [python, ./scripts/coverage.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          db_path:
            type: string
            format: project-path
          run_id:
            type: string
          collection_id:
            type: string
          md:
            type: boolean
          md_out:
            type: string
            format: project-path
        required: []
      approval: read
      output: ./schemas/coverage-result.json
      timeoutSeconds: 60
    - name: review_queue
      description: Diggable identity review queue (JSON; optional --md / --md_out).
      run: [python, ./scripts/review_queue.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          db_path:
            type: string
            format: project-path
          collection_id:
            type: string
          limit:
            type: integer
          no_shells:
            type: boolean
          no_versions:
            type: boolean
          md:
            type: boolean
          md_out:
            type: string
            format: project-path
          dig:
            type: boolean
          write_report:
            type: string
            format: project-path
        required: []
      approval: read
      output: ./schemas/review-queue-result.json
      timeoutSeconds: 60
    - name: versions
      description: detect | apply | propose | confirm version links.
      run: [python, ./scripts/versions.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          db_path:
            type: string
            format: project-path
          action:
            type: string
            description: detect|apply|propose|confirm
          collection_id:
            type: string
          from_id:
            type: string
          to_id:
            type: string
          relation:
            type: string
          confidence:
            type: number
          link_id:
            type: integer
          no_merge:
            type: boolean
          auto_merge_high_confidence:
            type: boolean
          evidence_json:
            type: string
          http:
            type: boolean
          providers:
            type: string
          mock:
            type: boolean
          fixture_dir:
            type: string
          no_wiki_repoint:
            type: boolean
          preview_fetch:
            type: boolean
        required: []
      approval: write
      output: ./schemas/versions-result.json
      timeoutSeconds: 180
---

# Literature library

One project, one `library.sqlite`. This skill finds, imports, and manages papers: collections, frozen search runs, HTTP or mock search, file ingest, deduplication, enrich, seed expansion, query, export, coverage, run close, identity review, and preprint–journal versions. Wiki prose is the `lit-wiki` skill.

Call the `lit_*` tools below. The old per-command skill folders (`lit-search`, `lit-ingest`, `lit-db`, and the rest of the library set) are this one skill. Script short names are unchanged, so the plugin prefix `lit` still yields `lit_search`, `lit_ingest`, `lit_collection`, and the other ids callers already use.

Do not use this skill for citation-verify, hypothesis cards, scientific-gate, or omics analysis. Do not add commercial gateways, Google Scholar scraping, embeddings, or a plugin API key.

## Tool order

1. `lit_collection` — create or list a 课题. Create initializes `library.sqlite` and, by default, `collections/<id>/wiki-vault/`. Set `no_wiki` only when the user does not want a vault.
2. `lit_run_config` — freeze `runs/<run_id>/run.yaml` and `run.json` (query, sources, dates, limit, provider notes, optional gold set and out-of-window controls). Open sources only, at most 8. `no_commercial_sources` stays true. The caller checks that gold DOIs are real; the plugin only checks the format.
3. `lit_search` — read that run, write `hits.jsonl`, and auto-ingest with dedup. `mock` is for offline smoke. `skip_ingest` keeps hits without writing the library. Optional CLI `--close_run` applies the same gate as `lit_run_close`.
4. `lit_ingest` — hits JSONL, BibTeX, RIS, a DOI/PMID list, a PDF text sidecar, or an authorized local PDF. ID-list shells auto-enrich unless `no_enrich`.
5. `lit_enrich` — fill title, abstract, authors, year, and venue on ID-only shells before any wiki compile. `mock` is offline. Richest-merge only; never invent a field.
6. `lit_expand` — similar, citations, or references from a seed. **Preview only** (`wrote: false`) until `confirm_write`. Do not call the graph closed unless the result says `saturated`.
7. `lit_query` — FTS5 on title, abstract, and tags, or `LIKE` when FTS is unavailable. No embeddings.
8. `lit_export` — `search_log.jsonl` for research-evidence. This plugin does not write hypothesis cards or run a gate.
9. `lit_coverage` and `lit_run_close` — report coverage, then mark the run `passed` or `needs_attention`. `fail_under_coverage` exits 2 when the gate fails.
10. `lit_review_queue` — identity adjudication only, not full-text review. `dig` is the schedulable form: an empty backlog exits 0 with `silent: true` and writes no report. Do not create a Grok Bot routine in this plugin.
11. `lit_versions` — preprint↔journal links. `preview_fetch` creates a metadata shell for a missing opposite DOI and does not pull a PDF. Confirm merges into the journal DOI; preprint ids stay aliases. Wiki pages are repointed unless `no_wiki_repoint`.

`lit_init_db` creates or upgrades the schema (idempotent). `lit_merge_ids` plans complementary DOI/PMID merges and writes them only when `apply` is set. Ingest already runs the high-confidence pass.

## Dedup

Strong ids merge in this order: DOI > PMID > PMCID > arXiv > OpenAlex. The publication primary key is that natural key (`doi:…`, `pmid:…`, …), not a UUID. A later higher-priority id reassigns the row; dependents and aliases follow.

A soft fingerprint (title + first author + year) or a real strong-id clash becomes `needs_review`. Do not auto-merge those. Complementary strong-id shells of the same work merge toward the DOI. An ids-file line (or consecutive lines) that carries both a DOI and a PMID is one publication.

## Search

Default backend is in-plugin HTTP (`scripts/providers_http.py`). This plugin does not ship an MCP server and does not bundle API keys.

| Source | Behavior |
| --- | --- |
| pubmed, europepmc, openalex, crossref | live public REST |
| arxiv | Atom; q-bio first, then `all:` |
| biorxiv, medrxiv | no keyword API: DOI detail, a date window, or Europe PMC `JOURNAL` |
| pmc | E-utilities `db=pmc` |
| doaj, datacite, zenodo | live public APIs; a shape change soft-fails |
| semantic_scholar | live only with host `S2_API_KEY` (never stored); 401/403/429 returns empty hits plus meta |
| mock | two fake hits, no network |

Allowed source names: pubmed, europepmc, openalex, crossref, semantic_scholar, arxiv, biorxiv, medrxiv, pmc, doaj, datacite, zenodo.

Optional `--backend mcp` or `LIT_SEARCH_BACKEND=mcp` calls a host `LIT_MCP_CMD` bridge, or `--mcp_mock_fixture` in tests. If MCP is down, soft-fail can fall back to HTTP. Do not add CORE, Scopus, Web of Science, SmartLib, CNKI, or Google Scholar.

Optional polite-pool mailto: `OPENALEX_MAILTO` or `CONTACT_EMAIL`.

## Files, OA, and campus access

Digital PDFs use `pdftotext` when it is installed. OCR (`tesseract` or `LIT_PDF_OCR_CMD`) is an optional soft-fail. TextIn is not required and is not the knowledge base.

Ingest records `pub_files` (path, `content_hash`, `ingest_origin`). The path is not the identity. `lit_rebind_files` matches `content_hash` after a move and stays a dry-run until `confirm`.

`scripts/watch_rebind.py` is the same rebind on a stdlib poll (`--once` or `--interval`). It is a script in this skill, not its own tool: `skill_run` it as `watch_rebind.py`, or run `python skills/literature-library/scripts/watch_rebind.py`. Dry-run unless `--confirm`.

`--verify_oa` classifies a candidate URL as `verified_oa`, `unknown`, `paywall_or_blocked`, or `not_pdf` and stores that on the publication. Never save paywalled bytes as library full text, and never treat an unverified OA link as an owned PDF.

`lit_institutional_access` requires `institutional_access: true`.

| method | What to do |
| --- | --- |
| `user_provided_file` | mount an authorized local PDF (preferred) |
| `vpn_proxy` | fetch through `LIT_INSTITUTIONAL_PROXY`, then mount |
| `session_reuse` | the environment already has the user's VPN or session |
| `saml_browser` | do not start a browser here; the host finishes SAML, then use `user_provided_file` |

Only the user's own campus account. Credentials stay in the user vault. Never write them into SQLite, and never use a borrowed login. Schema v6 stores `access_origin` and `access_method` on `pub_files`.

## Expand

| relation | First path | Fallback |
| --- | --- | --- |
| similar | PMID → PubMed elink `pubmed_pubmed` | OpenAlex `related_works` |
| references | OpenAlex `referenced_works` | OpenCitations, then OpenAlex by DOI |
| citations | OpenAlex `cites:` | OpenCitations, then OpenAlex by DOI |

`max_rounds` defaults to 1. Two rounds in a row under `low_yield_ratio` (default 0.05) stop as saturated. Stopping because the round cap was hit is `max_rounds_not_saturated`, not closure. `confirm_write` writes `collections/<id>/expands/<exp-id>/hits.jsonl` and ingests. `dry_run` forces a preview even when confirm is set.

## Review and versions

Review reason codes:

| code | Meaning |
| --- | --- |
| `soft_fingerprint_conflict` | same title + first author + year, different strong ids |
| `conflicting_strong_ids` | two existing pubs or a `#conflict-` row |
| `empty_unenriched_shell` | id shell with no real title or abstract |
| `pending_version_link` | `publication_versions.needs_confirm=1` |

Version relations are `preprint_of`, `published_as`, and `same_work`. Detect is offline heuristics plus optional OpenAlex/Crossref HTTP (`http`, or `mock` fixtures). `apply` writes pending links; `auto_merge_high_confidence` merges into the journal DOI. There is no SmartLib or other commercial version source.

Coverage also reports gold-set recall (`evidence_weak` when the gold set is empty or shorter than 5), out-of-window leakage, and declared sources that were not contacted. Hints are honest gaps, not a made-up recall number.

## Shared modules

`scripts/db.py` (schema v6), `oa_verify.py`, `providers_http.py`, `providers_mcp.py`, and `pdf_extract.py` are libraries the CLIs import. Use the tools above rather than reimplementing dedup or HTTP.

## Deliver

Reply with absolute paths and the inserted / merged / `needs_review` counts. Defaults live in the project: `library.sqlite`, `collections/<id>/runs/`, `collections/<id>/wiki-vault/`, `exports/search_log.jsonl`.
