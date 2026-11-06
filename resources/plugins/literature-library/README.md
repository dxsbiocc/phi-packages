# literature-library（文献库）· 0.6.1

Project-scoped literature library in `phi-packages` at `resources/plugins/literature-library/`. **0.6.1** groups the tools into two skills, `literature-library` and `lit-wiki`. Script CLIs and `lit_*` tool ids are unchanged from the 0.6.0 landing. Live `demo-*` trees are not part of this package; offline tests live under `tests/`.

English: Project-scoped literature library for Phi — HTTP search by default, optional **phi-literature MCP client** (this plugin does not ship an MCP server), file ingest, dedup, SQLite identity store + **LLM Wiki** compiled knowledge layer. **No embeddings. No commercial literature gateways. No plugin API keys.** Does **not** do citation-verify, hypothesis cards, or scientific-gate. TextIn is not the KB backbone.


## 最短真实工作流 / Shortest real workflow

```bash
# 1) collection (+ wiki-vault)
python skills/literature-library/scripts/collection.py --project_root P --name "topic" --create

# 2) freeze run (+ optional gold_set)
python skills/literature-library/scripts/run_config.py --project_root P --collection_id C \
  --query "…" --sources pubmed,openalex --gold_set tests/fixtures/gold.ids

# 3) search → enrich → expand → coverage/close
python skills/literature-library/scripts/search.py --project_root P --collection_id C --run_id R --limit 20
python skills/literature-library/scripts/enrich.py --project_root P --collection_id C
python skills/literature-library/scripts/expand.py --project_root P --collection_id C --seed_ids … --max_rounds 2 --confirm_write
python skills/literature-library/scripts/coverage.py --project_root P --run_id R --md
python skills/literature-library/scripts/run_close.py --project_root P --run_id R

# 4) wiki: agent compile (NOT mock for live)
python skills/lit-wiki/scripts/wiki_agent_status.py --project_root P --collection_id C
python skills/lit-wiki/scripts/wiki_compile_request.py --project_root P --collection_id C --publication_id ID --ensure_vault
# → you fill compile_response.json with your model
python skills/lit-wiki/scripts/wiki_compile_apply.py --project_root P --collection_id C --publication_id ID --response compile_response.json
python skills/lit-wiki/scripts/concept_merge_request.py --project_root P --collection_id C --slug t-cell-exhaustion
# → fill merge response → concept_merge_apply
python skills/lit-wiki/scripts/wiki_lint.py --project_root P --collection_id C
python skills/lit-wiki/scripts/wiki_query.py --project_root P --collection_id C --query "exhaustion"

# 5) optional: local files moved → rebind by hash (or path-watch)
python skills/literature-library/scripts/rebind_files.py --project_root P --scan_dir ./papers   # dry-run
python skills/literature-library/scripts/watch_rebind.py --project_root P --scan_dir ./papers --once --confirm

# 6) optional: institutional authorized PDF (own account only)
python skills/literature-library/scripts/institutional_access.py \
  --project_root P --collection_id C --institutional_access \
  --method user_provided_file --file_path ./papers/authorized.pdf --doi 10.…
```

**Wiki default = agent compile.** `--compile mock` is for offline CI only — do **not** recommend mock for live acceptance.

## Architecture / 架构

```
User / LiteratureLibrary agent
        │
        ├─ lit-collection  →  collections + library.sqlite + wiki-vault/ (default)
        ├─ lit-run-config  →  runs/<id>/run.yaml + frozen run.json
        ├─ lit-search      →  hits.jsonl (+ hits.meta.json) → auto lit-ingest
        ├─ lit-run-close   →  coverage.json + gate → passed|needs_attention
        │       ├─ draft: in-plugin HTTP (pubmed/…/medrxiv/pmc/doaj/datacite/zenodo; S2 soft-fail) / --mock
        │       └─ optional: --backend mcp / LIT_SEARCH_BACKEND=mcp → providers_mcp client (plugin does NOT ship MCP server)
        ├─ lit-ingest      →  bib / RIS / IDs / PDF (pdftotext + optional OCR) → dedup upsert
        ├─ lit-institutional → own-account campus/SSO authorized fulltext mount (explicit flag)
        ├─ lit-enrich      →  ID-shell → PubMed/EuropePMC title+abstract (auto after ids)
        ├─ lit-query       →  FTS5 / LIKE (SQLite)
        ├─ lit-export      →  search_log.jsonl → research-evidence
        ├─ lit-expand      →  HTTP related/citations/references (confirm_write)
        ├─ lit-review-queue → identity human_queue (soft/conflict/shell/version)
        ├─ lit-coverage    → run/collection coverage report
        ├─ lit-versions    → preprint↔journal links (confirm / merge)
        └─ Wiki layer (LLM Wiki; no embeddings)
                ├─ lit-wiki-init / ingest / lint / query
                ├─ concept_merge_request → agent JSON → concept_merge_apply
                └─ lit-wiki-suggest-tree → apply-tree --confirm
```

Vault per collection:

```
collections/<id>/wiki-vault/
  AGENTS.md          # schema
  raw/publications/  # immutable snapshots
  wiki/
    index.md         # theme tree
    log.md
    sources/ concepts/ entities/
```

## Boundaries / 边界

| | literature-library | phi-literature MCP | research-evidence | paper-repro |
|--|--|--|--|--|
| Role | Store + dedup + SQLite + Wiki | Online multi-source search | Evidence / hypotheses / gate | Deposit reproduction |
| Owns DB / vault | Yes | No | Consumes export | May read |
| MCP path | — | `/workspace/phi-packages/mcp-servers/literature` | — | — |

**phi-literature** is a first-party Phi MCP package — **not** a Claude plugin.  
**TextIn** document-tree is **not** the KB backbone (optional PDF parser later; not integrated).

## Skills（0.6.1，toolPrefix `lit`）

Two skills. Scripts are the tools; each skill teaches one capability. Every former `lit_*` tool id is unchanged.

### `literature-library` — find / import / manage

| Tool | Role |
|--|--|
| `lit_init_db`, `lit_merge_ids` | schema and complementary strong-ID repair (`scripts/db.py`) |
| `lit_collection` | create or list a 课题; default wiki-init |
| `lit_run_config` | freeze `run.json` |
| `lit_search` | in-plugin HTTP (pubmed through zenodo; S2 soft-fail) plus mock; optional MCP client; optional `--close_run` |
| `lit_ingest`, `lit_rebind_files` | files, `pdf_extract.py`, auto-enrich unless `--no_enrich`. `watch_rebind.py` is a script in this skill, not its own tool |
| `lit_institutional_access` | own-account full text; requires `--institutional_access` |
| `lit_enrich` | HTTP plus `--mock` |
| `lit_expand` | HTTP plus `--mock`; write only with `--confirm_write` |
| `lit_query` | FTS or LIKE |
| `lit_export` | `search_log.jsonl` |
| `lit_run_close` | coverage gate and run status |
| `lit_coverage` | run or collection coverage |
| `lit_review_queue` | `--dig`; suggested actions |
| `lit_versions` | detect / apply / propose / confirm; `--preview_fetch` |

### `lit-wiki` — LLM Wiki compile

| Tool | Role |
|--|--|
| `lit_wiki_init` | vault skeleton |
| `lit_wiki_ingest` | agent compile entry |
| `lit_wiki_compile_request`, `lit_wiki_compile_apply` | request, then the model’s JSON, then write |
| `lit_concept_merge_request`, `lit_concept_merge_apply` | cross-paper concept merge |
| `lit_wiki_agent_status` | missing sources and pending merges |
| `lit_wiki_query` | keyword search; no vectors |
| `lit_wiki_lint` | vault health |
| `lit_wiki_suggest_tree` | proposals only |
| `lit_wiki_apply_tree` | write only with `--confirm` |

`wiki_lib.py` and `wiki_llm.py` are helpers those tools import. They are not separate tools.

## Why heuristic existed (v0.2)

v0.2 used deterministic title/abstract heuristics so smoke tests and offline demos worked without calling a model. That scaffolding is **not** target LLM Wiki quality.

## Wiki compile (v0.4 — invoking agent's model)

**Normal use needs no API key.** Wiki prose is produced by the **LiteratureLibrary / host agent's already-running model**. Plugin scripts only build grounded prompts + JSON schema, validate the agent's JSON, and write the vault.

```bash
# 1) request (default --compile agent)
python skills/lit-wiki/scripts/wiki_compile_request.py \
  --project_root P --collection_id C --publication_id ID --ensure_vault

# 2) agent fills compile_response.json with its own model

# 3) apply
python skills/lit-wiki/scripts/wiki_compile_apply.py \
  --project_root P --collection_id C --publication_id ID \
  --response path/to/compile_response.json
```

- Offline smoke: `--compile mock`
- Emergency offline: `--compile heuristic` or `LIT_WIKI_COMPILE=heuristic`
- Do **not** configure plugin `OPENAI_API_KEY` / `LIT_WIKI_LLM_API_KEY` — that path was removed.

## Dedup / Publication PK (v0.4.6) + enrich (v0.4.2)

- **Publication PK** = natural key on `publications.id` (not UUID): `doi:` > `pmid:` > `pmcid:` > `arxiv:` > `openalex:` > `soft:` > `tmp:`. Schema v2 UNIQUE indexes on strong ID columns.
- **Canonical upgrade (0.4.6):** pmid-only (or pmcid/arxiv) row that later gains a higher-priority ID is **reassigned** to the preferred PK (dependents + aliases follow); soft fingerprint never auto-merges.
- **Ids-file co-identity (0.4.1):** consecutive/same-line `doi`+`pmid` in an `.ids` file → one pub with both aliases (no bib required).
- Priority: DOI > PMID > PMCID > arXiv > OpenAlex; soft fingerprint → `needs_review` (no auto-merge).
- **Complementary strong IDs** (doi-only shell + pmid-only shell + joint hit, or `#conflict-` cluster with compatible titles) auto-merge into DOI-preferred canonical; true ID clashes stay `needs_review`.
- Repair CLI (dry-run default): `python skills/literature-library/scripts/merge_ids_cli.py --project_root <proj>`; `--apply` writes.
- Ingest post-pass also runs high-confidence merges.
- **Auto-enrich shells (0.4.2):** after ids ingest, `lit-enrich` fills title/abstract/authors/year/venue via PubMed efetch / Europe PMC (richest-merge; never invents). Opt out: `--no_enrich`. Offline: `--enrich_mock` / `enrich.py --mock`.
- Prefer enrich before wiki `compile_request` when abstract missing (request JSON warns).

## Smoke

From this plugin root:

```bash
bash tests/run_smoke.sh
bash tests/run_smoke_v2.sh
```

## Environment

`phi:python@1` only. No private conda. No embedding deps.

## Search providers (0.4.3+ / 0.4.9)

The search CLI uses `skills/literature-library/scripts/providers_http.py` (public REST only; no keys in package).
Optional env: `OPENALEX_MAILTO` / `CONTACT_EMAIL`, `S2_API_KEY` (never committed).
`--sources openalex,crossref,arxiv` or `--sources pmc,doaj,datacite,zenodo` overrides `run.json`.
Live `demo-*` trees (for example `demo-search-v043/`) are not included in this package.

| source | endpoint family |
|--|--|
| pubmed / europepmc | NCBI / Europe PMC |
| openalex / crossref | OpenAlex / Crossref |
| arxiv / biorxiv / medrxiv | arXiv Atom / bioRxiv API (+ EPMC fallback) |
| **pmc** (0.4.9) | NCBI E-utilities `db=pmc` |
| **doaj** (0.4.9) | DOAJ `/api/search/articles` |
| **datacite** (0.4.9) | DataCite `/dois?query=` |
| **zenodo** (0.4.9) | Zenodo `/api/records` |
| semantic_scholar | Graph API (soft-fail without key) |

No CORE / Scopus / WoS / SmartLib / CNKI.

## Expand (0.4.4)

```bash
python skills/literature-library/scripts/expand.py \
  --project_root P --collection_id C \
  --seed_id 10.1016/j.cell.2017.05.035 \
  --relation similar --limit 5
# preview only → wrote:false

python skills/literature-library/scripts/expand.py ... --confirm_write   # ingest
```

| relation | primary | fallback |
|--|--|--|
| similar | PubMed elink (PMID) | OpenAlex `related_works` |
| references | OpenAlex `referenced_works` | OpenCitations → OpenAlex DOI hydrate |
| citations | OpenAlex `cites:W…` | OpenCitations → OpenAlex DOI hydrate |

Live `demo-expand-v044/` is not included in this package. No in-plugin MCP server; no new pip deps.

## Concept deep merge (0.4.5)

Per-source contributions: `wiki/concepts/_contrib/<slug>/<safe_id>.json` (deterministic).
Agent path: `concept_merge_request` → fill JSON → `concept_merge_apply` / `--llm-json`.
Offline: `--mock`. Ingest sets `concept_merge_pending` when a concept already had another source.
`needs_review` pubs require a banner on the concept page. Live `demo-concept-merge-v045/` is not included in this package.

## Human queue + coverage + version chain (0.4.7)

**No third-party commercial literature gateways** (no SmartLib / global-biblio-base). Public HTTP adapters only.

```bash
# Identity review queue (soft conflicts, ID clashes, shells, pending version links)
python skills/literature-library/scripts/review_queue.py \
  --project_root P --collection_id C --md

# Coverage for a search run / collection
python skills/literature-library/scripts/coverage.py \
  --project_root P --run_id run-... --md

# Preprint↔journal version links
python skills/literature-library/scripts/versions.py --project_root P --action detect
python skills/literature-library/scripts/versions.py --project_root P --action apply
python skills/literature-library/scripts/versions.py --project_root P --action apply \
  --auto_merge_high_confidence   # merge into journal DOI; preprint IDs stay as aliases
```

Schema v3: `publication_versions` + optional `publications.review_reasons_json`.

## HTTP version detect + wiki repoint (0.4.8)

```bash
# Offline fixtures (no network)
python skills/literature-library/scripts/versions.py --project_root P --action detect --mock

# Live OpenAlex / Crossref (soft-fail; public HTTP only)
python skills/literature-library/scripts/versions.py --project_root P --action apply \
  --http --auto_merge_high_confidence

# After id upgrade/merge, wiki/sources/<old_safe>.md → <new_safe>.md (+ stub)
# Automatic on lit-versions apply/confirm and lit-db merge_ids --apply
# (pass --no_wiki_repoint to skip)
```

Schema unchanged (still v3). No commercial gateways; no MCP; no new deps.

## Remaining gaps (honest)

- Soft title-fingerprint conflicts still need human review (by design); use `lit-review-queue`.
- Entity extraction is lexicon/heuristic — rare genes outside the allowlist may be missed; some allcaps acronyms may still slip through.
- OpenAlex/Crossref version ties auto-merge when both sides exist; missing opposite DOI can be **preview-fetched** (metadata-only shell + pending link; no PDF) then confirmed.
- Institutional full text is an own-account scaffold (`lit-institutional`); no CNKI/Scholar scrape. Live SSO/proxy soak is still open. File-ingest OA verify is soft classify only (0.5.1).
- Live HTTP version detect depends on provider uptime/rate limits (soft-fail; use `--mock` offline).
- Concept `_contrib` / index links are string-repointed; deep concept page prose may still mention old ids until re-merge.
- Optional MCP *client* (`--backend mcp`) needs a host `LIT_MCP_CMD` (or a mock fixture). HTTP stays the default. This plugin does not ship an MCP server.
- Agent-quality concept merges still depend on the host model filling JSON (scripts only validate/apply; mock is for tests).


## Open search sources (0.4.9)

Added free HTTP providers: **pmc**, **doaj**, **datacite**, **zenodo**.

```bash
python skills/literature-library/scripts/search.py \
  --project_root P --collection_id C --run_id R \
  --sources pmc,doaj,datacite,zenodo
```


## Protocol freeze + ops (0.5.0)

```bash
# 1) Freeze protocol
python skills/literature-library/scripts/run_config.py \
  --project_root P --collection_id C --query "…" \
  --sources pubmed,openalex,pmc --date_from 2018-01-01 --limit 20

# 2) Search (+ optional gate)
python skills/literature-library/scripts/search.py \
  --project_root P --collection_id C --run_id R --mock \
  --close_run --min_unique_pubs 5

# Or close separately
python skills/literature-library/scripts/run_close.py \
  --project_root P --run_id R --min_unique_pubs 5 --fail_under_coverage

# 3) Schedulable human dig (silent if empty)
python skills/literature-library/scripts/review_queue.py \
  --project_root P --collection_id C --dig --write_report /tmp/lit-dig.md

# 4) Version chain: preview-fetch missing opposite DOI (metadata only)
python skills/literature-library/scripts/versions.py \
  --project_root P --action apply --mock --preview_fetch
python skills/literature-library/scripts/versions.py \
  --project_root P --action confirm --link_id 1
```

Schema v4: `runs.status` / `coverage_path` / `gate_json`. No commercial sources; no MCP; no new pip deps.


## File ingest polish (0.5.1)

```bash
# PDF text sidecar → pub_files + optional OA verify (mock offline)
python skills/literature-library/scripts/ingest.py \
  --project_root P --collection_id C \
  --input_path tests/fixtures/pdf_sidecar_zheng.txt --format pdf_text \
  --verify_oa --oa_mock_fixture tests/fixtures/oa_mock_head.json

# Authorized local PDF (sidecar or filename IDs only; no binary parse)
python skills/literature-library/scripts/ingest.py \
  --project_root P --collection_id C \
  --input_path ./papers/10.1016_j.cell.2017.05.035.pdf --format pdf

# Enrich + OA verify
python skills/literature-library/scripts/enrich.py \
  --project_root P --collection_id C --mock --verify_oa \
  --oa_mock_fixture tests/fixtures/oa_mock_head.json
```

Schema **v5**: `pub_files` + `publications.oa_status`/`oa_json`. Candidate OA URLs are never stored as owned PDF bytes.


## Recall evidence + snowball saturation (0.5.2)

```bash
# Freeze gold-set + out-of-window controls with the run protocol
python skills/literature-library/scripts/run_config.py \
  --project_root P --collection_id C --query "…" \
  --sources pubmed,openalex,doaj --date_from 2018-01-01 --date_to 2024-12-31 \
  --gold_set "10.1000/a,10.1000/b,pmid:123,10.1000/c,10.1000/d" \
  --out_of_window_controls "10.1000/old.paper"

# Coverage includes gold recall / OOW leakage / uncontacted sources
python skills/literature-library/scripts/coverage.py \
  --project_root P --run_id R --md

# Snowball with honest saturation (never false "closed")
python skills/literature-library/scripts/expand.py \
  --project_root P --collection_id C --seed_id doi:10.1000/seed \
  --relation similar --max_rounds 5 --mock
# → saturated + stop_reason, or max_rounds_not_saturated (closure NOT claimed)
```

Schema unchanged (v5). Inspiration: research-skills lit-search recall/saturation honesty (conceptual only).



## 0.6.0 highlights

- **PDF extract:** `pdftotext` when available; optional OCR soft-fail (`tesseract` / `LIT_PDF_OCR_CMD`). TextIn is not the KB backbone.
- **Path watch:** `watch_rebind.py` polls and rebinds by `content_hash`.
- **Institutional access:** `lit-institutional` (explicit `--institutional_access`; own VPN/proxy/file only).
- **MCP client (optional):** `--backend mcp` + `LIT_MCP_CMD` / `--mcp_mock_fixture`; HTTP remains default. No MCP server is shipped in this plugin.
- **research-evidence:** optional `--delegate_lit` / `--hits_jsonl` bridge (plugins stay separate; not merged here).
- **Schema v6:** `pub_files.access_origin` / `access_method`.
- See `NOTES-0.6.0.md` + `PACKAGING.md`. This tree is the first packaged landing (`version: 0.6.0`).
