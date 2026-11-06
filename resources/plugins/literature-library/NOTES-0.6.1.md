# Notes — literature-library 0.6.1

Packaging only. Script CLIs, SQLite schema, wiki contracts, and `lit_*` tool ids are unchanged from 0.6.0. No commercial gateways, Scholar scrape, embeddings, or plugin API keys were added.

Version is a patch: callers that invoke `lit_search`, `lit_ingest`, `lit_wiki_query`, and the other prefixed tools keep the same names. What changed is which skill directory teaches those tools.

## Why two skills

Phi skills should describe a capability. The 0.6.0 tree had one skill per CLI, which made the LiteratureLibrary agent load 21 skills to do two jobs. Sibling plugin `visualization` already groups catalog / profile / prepare / render under `omics-visualization` while `toolPrefix: viz` keeps tool ids stable. The same split here:

1. **`literature-library`** — find, import, and manage the project library.
2. **`lit-wiki`** — LLM Wiki compile layer (the invoking model writes prose; scripts validate and write the vault).

Review queue, coverage, run close, versions, enrich, expand, and institutional access stay in (1). They are steps in managing the library, not a third product surface. A third skill would only re-split that workflow.

## Old skill → new skill

| Former skill | Now |
| --- | --- |
| `lit-db`, `lit-collection`, `lit-run-config`, `lit-search`, `lit-ingest`, `lit-institutional`, `lit-enrich`, `lit-expand`, `lit-query`, `lit-export`, `lit-run-close`, `lit-coverage`, `lit-review-queue`, `lit-versions` | `skills/literature-library` |
| `lit-wiki-lib`, `lit-wiki-init`, `lit-wiki-ingest`, `lit-wiki-lint`, `lit-wiki-query`, `lit-wiki-suggest-tree`, `lit-wiki-apply-tree` | `skills/lit-wiki` |

`db.py` stays beside the library CLIs. Wiki CLIs import it from `skills/literature-library/scripts/`. Library CLIs that repoint wiki pages import `skills/lit-wiki/scripts/wiki_lib.py`. `watch_rebind.py` remains a script on the library skill and is still not its own tool.

## Counts

`counts.pluginSkills` is 4 (visualization 2 + literature-library 2).
