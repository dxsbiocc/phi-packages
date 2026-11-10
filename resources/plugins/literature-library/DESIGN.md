# Design — literature-library 0.7.0

## Scope

A project keeps its literature in one place: **find it, store it, don't duplicate it, search it, cite it,
and write it up.** Anything else belongs elsewhere.

| In | Out |
|--|--|
| Search open sources and save hits with deduplication | The search services themselves (a separate literature MCP server and the PubMed connector exist) |
| Import `.bib`, `.ris`, ID lists, PDFs, folders | Fetching paywalled full text, handling credentials, institutional access brokering |
| Preprint <-> journal links, identity review, moved-file repair | Citation fact-checking, evidence scoring, hypotheses, omics analysis |
| Export BibTeX, RIS, CSL-JSON, CSV | Reference-style formatting (citeproc does that) |
| A wiki the agent writes from the saved text | Embeddings, keyword-list "NLP", plugin-side LLM calls |

## Decisions

1. **Skills are workflows, tools are actions.** Two skills, 12 tools. The agent file holds only role and
   boundaries. (The previous 21 skills mirrored 27 scripts one to one.)
2. **One skill owns the code.** `manage-library` carries every script; `compile-wiki` is documentation
   that points at its tools. The plugin contract has no shared-code directory, and importing between
   skills made none of them installable alone.
3. **SQLite for identity, markdown for reading.** `library.sqlite` decides which records are one paper.
   The wiki vault never decides identity; merging ids renames wiki pages to follow.
4. **External vs internal search are different tools.** `lit_search` reaches the network; `lit_query`
   never does. Both descriptions say so, because models confuse them otherwise.
5. **Batch the agent's writing.** Wiki prose comes from the invoking agent's own model. One request file
   covers up to 25 papers and one apply call writes them; each result is validated separately so one
   malformed entry does not lose the rest.
6. **No hand-rolled language processing.** The old keyword lists for methods, diseases and gene symbols
   were tuned to one oncology demo and needed noise blocklists. The agent reads the abstract; the script
   checks the shape of what it wrote.
7. **Preview before writing outward-facing changes.** `lit_expand` previews unless `confirm_write`;
   `merge_ids` and `rebind` are dry runs unless `apply` / `confirm`.
8. **Uncertain identity stays uncertain.** Title + first author + year alone never merges two records;
   they become `needs_review` and wait for a person.

## Search record

A search is a `run`: query, sources, dates, per-source counts and time, in
`collections/<id>/runs/<run_id>/run.json`, with the raw `hits.jsonl` beside it. That is enough to
describe a search in a methods section. Pass/fail gates, gold-set recall and saturation claims were
removed: they answered systematic-review questions this library does not own.

## Schema (SQLite, version 6)

`collections`, `publications`, `aliases`, `collection_pubs`, `runs`, `run_hits`, `pub_files`
(path, content hash, ingest origin, access origin and method), `publication_versions`, `meta`, and an
FTS5 table `pubs_fts` (LIKE fallback when FTS5 is missing). `theme_nodes` and `theme_pubs` remain in the
schema for old databases but are no longer written: themes live in `wiki/index.md`.

## Sources

`pubmed`, `europepmc`, `openalex`, `crossref`, `arxiv`, `biorxiv`, `medrxiv`, over public REST with the
standard library. Sources that needed keys or served other kinds of records (Semantic Scholar, PMC,
DOAJ, DataCite, Zenodo) were dropped; add one back only with a concrete need.

## Known limits

- Live HTTP paths are covered by fixtures, not by CI against the real services.
- A PDF with no identifier (a thesis, a proposal, a report) is identified by its content hash, `tmp:<sha256 prefix>`: the same file imported again, or under another name, is the same record. Two differently typeset copies of one paper are not detected as duplicates unless a registry id or title match joins them.
- Wiki quality depends on the abstract; `with_fulltext` adds up to the first N characters of a saved PDF.
