# Changelog

## 0.7.0

Restructured around tasks instead of scripts. **Tool names changed; there are no aliases.**

### Layout

- 21 skills -> 2: `manage-library` (all scripts and tools) and `compile-wiki` (workflow only).
- 27 script tools -> 12. Internal modules moved to `skills/manage-library/scripts/lib/`; no skill
  imports another skill's code.
- `LiteratureLibrary` agent shortened to a role and boundaries; workflows live in the skills.

### Tool mapping

| Before | Now |
|--|--|
| `lit_run_config` + `lit_search` + `lit_run_close` | `lit_search` (takes `query` and `sources`, creates the run) |
| `lit_ingest`, `lit_institutional` | `lit_import` (`access_origin` records how a file was obtained) |
| `lit_enrich`, `lit_review_queue`, `lit_versions`, `lit_db` merge, `lit_rebind_files` | `lit_curate` (`action`) |
| `lit_query`, `lit_wiki_query` | `lit_query` (wiki hits included) |
| `lit_export` (search_log only) | `lit_export` (`bibtex`, `ris`, `csl-json`, `csv`, `search_log`) |
| `lit_wiki_ingest`, `lit_wiki_compile_request/apply`, `lit_wiki_agent_status` | `lit_wiki_pending` + `lit_wiki_apply` (one batch of up to 25 papers) |
| `lit_wiki_init` | done by `lit_collection` and `lit_wiki_pending` |
| `lit_wiki_lint` | `lit_wiki_lint` |
| `lit_coverage`, `lit_concept_merge_*`, `lit_wiki_suggest_tree`, `lit_wiki_apply_tree` | see "Removed" |

### Added

- `lit_status`: one-call overview with `next_steps`.
- `lit_get`: full record of one paper, optional PDF text and wiki page.
- `lit_export` formats BibTeX, RIS, CSL-JSON, CSV with stable cite keys.
- `lit_import` accepts a folder; `lit_collection list` shows paper counts.
- `collection_id` may be omitted when the project has one collection.
- Wiki pages can be written from saved PDF text (`lit_wiki_pending with_fulltext`).

### Removed

- Coverage gates, `gold_set`, out-of-window controls, saturation criteria, scheduled dig mode.
- Institutional access scaffold, path-watch daemon, MCP client backend.
- Sources `pmc`, `doaj`, `datacite`, `zenodo`, `semantic_scholar`.
- Keyword-list extraction of methods, diseases and genes; offline `mock` and `heuristic` compile
  modes; the concept-merge request/apply pipeline; theme-tree suggest/apply (edit `index.md`).
- Entity pages: genes and proteins are now ordinary concepts. Existing `entities/` pages are still
  linted.

### Compatibility

`library.sqlite` schema is unchanged (v6) and old databases open as before. Old wiki vaults still
lint; `wiki/concepts/_contrib/` files are no longer read.
