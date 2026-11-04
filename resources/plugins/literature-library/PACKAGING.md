# Packaging — literature-library 0.6.0

First package landing in phi-packages:

```
resources/plugins/literature-library/
```

Version is `0.6.0` (`-draft` dropped). This file records what was included. It is not a block on further edits.

## Included

| Path | Notes |
|--|--|
| `phi-package.yaml` | `id: literature-library`, `type: plugin`, `toolPrefix: lit`, `version: 0.6.0` |
| `agents/` | LiteratureLibrary.md |
| `skills/` | all `lit-*` skills, including `lit-institutional` |
| `tests/` | offline suite and fixtures |
| `DESIGN.md` `README.md` `REVIEW-CHECKLIST.md` `NOTES-0.6.0.md` `PACKAGING.md` `UX-NOTES-0.5.3.md` | docs |

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
2. `version: 0.6.0`
3. `components.skills` lists every skill directory
4. Summary stays within the plugin contract (300 characters)
5. Host `validatePlugin` is not in this repository; the contract check is local. Offline tests are under `tests/`
