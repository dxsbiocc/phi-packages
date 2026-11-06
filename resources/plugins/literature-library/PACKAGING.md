# Packaging — literature-library 0.6.1

```
resources/plugins/literature-library/
```

0.6.0 was the first landing in phi-packages (`-draft` dropped). **0.6.1** only regroups skills: 21 one-script directories become two capability skills. Tool ids stay `lit_*`. This file records what is included. It is not a block on further edits.

## Included

| Path | Notes |
|--|--|
| `phi-package.yaml` | `id: literature-library`, `type: plugin`, `toolPrefix: lit`, `version: 0.6.1` |
| `agents/` | LiteratureLibrary.md loads `literature-library` and `lit-wiki` |
| `skills/literature-library` | library lifecycle: db, collection, run config, search, ingest, institutional access, enrich, expand, query, export, coverage, run close, review queue, versions |
| `skills/lit-wiki` | wiki compile: vault helpers, init, ingest, compile request/apply, concept merge, lint, query, suggest-tree, apply-tree |
| `tests/` | offline suite and fixtures |
| `DESIGN.md` `README.md` `REVIEW-CHECKLIST.md` `NOTES-0.6.0.md` `NOTES-0.6.1.md` `PACKAGING.md` `UX-NOTES-0.5.3.md` | docs |

Script `name` fields are the same short names as 0.6.0, so prefixed tool ids do not change. See `NOTES-0.6.1.md` for the old-skill → new-skill map.

## Kept out

| Path | Why |
|--|--|
| `demo-*` / `demo-acceptance-*` / `demo-project*` source trees | bulky live demos; not part of this package |
| `.git` / `__pycache__` / `*.pyc` | noise |
| Vault credentials / personal PDF corpora | never ship |
| `research-evidence` | separate plugin; do not merge |
| Commercial gateways, CNKI, Google Scholar scrape, paywall bypass | product constraint |
| Embeddings, plugin API keys, an in-plugin MCP server | product constraint |

`tests/run_smoke.sh` may create a temporary `demo-project/` under the plugin root while it runs. That directory is a run artifact, not shipped source.

## Checks at landing

1. `toolPrefix: lit` does not collide with `viz` or standalone skill prefixes in this repo
2. `version: 0.6.1`
3. `components.skills` lists both skill directories and no others
4. Summary stays within the plugin contract (300 characters)
5. Host `validatePlugin` is not in this repository; the contract check is local. Offline tests are under `tests/`
