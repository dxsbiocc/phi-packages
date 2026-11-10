# Searching external sources

`lit_search` queries each chosen source once and saves what comes back. Each source has
its own query language and limits.

| Source | Good for | Query notes |
|--|--|--|
| `pubmed` | Biomedical literature | MeSH and field tags work: `"T-Lymphocytes"[MeSH] AND exhaustion`, `[ti]`, `[au]`. Dates filter by year only. Newest first. |
| `europepmc` | Biomedical plus preprints and some full-text records | Plain words or Europe PMC syntax. Full dates work. |
| `openalex` | Broad coverage, any field, citation links | Plain words. Full dates work. |
| `crossref` | DOI-registered works, books, conference papers | Plain words. Ranked by relevance, not recency. |
| `arxiv` | Physics, maths, computer science, q-bio | Dates are ignored by this source. |
| `biorxiv`, `medrxiv` | Preprints | No keyword API: a DOI returns that record, dates browse a window, and a plain topic goes through Europe PMC's journal filter. |

## Choosing sources

- Default is `pubmed,europepmc`. They overlap heavily; the library merges duplicates.
- Add `openalex` when recall matters or the topic is not purely biomedical.
- Add `biorxiv` or `arxiv` only when the user wants preprints.
- At most four sources per call. If you need more, make a second call.

## Reading the result

- `per_source[i].ok = false` means that source failed (network, rate limit, odd response). The
  other sources still ran and their hits are saved. Say which one failed and offer to retry it.
- `ingest.counts`: `inserted` are new papers, `merged` matched a paper already saved,
  `needs_review` have an uncertain identity and wait for `lit_curate action=queue`.
- The raw hits of every search are kept in `collections/<id>/runs/<run_id>/hits.jsonl` next to
  `run.json`, which records the query, sources, dates and when it ran. Use them as the search
  record when writing a methods section.

## Politeness

OpenAlex and Crossref ask for a contact address. If the user sets `OPENALEX_MAILTO`,
`CROSSREF_MAILTO` or `CONTACT_EMAIL` in the environment, it is sent; nothing is packaged.
Keep `limit` at the default unless the user asks for more.
