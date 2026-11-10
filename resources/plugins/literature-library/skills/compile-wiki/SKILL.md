---
name: compile-wiki
description: >-
  Turn saved papers into a readable wiki: write each paper's summary, key
  claims and concepts from its abstract (or saved PDF text), write concept
  pages that span several papers, and group the index into themes. Use after
  papers are in the library, for 整理文献, 文献笔记, 生成文献 wiki, 概念页,
  主题整理, 综述提纲. Not for searching, importing or exporting papers; that is
  manage-library.
phi:
  attachTo:
    - LiteratureLibrary
---

# Compile Wiki

Each collection has a vault at `collections/<id>/wiki-vault/` beside `library.sqlite`. This
skill fills it. **You write the prose with your own model**; the tools only prepare the
request, check what you wrote, and put the files in place. The tools belong to
`skill://manage-library`: `lit_wiki_pending`, `lit_wiki_apply`, `lit_wiki_lint`.

Page formats and the vault layout: `references/pages.md`.

## Write source pages (one batch at a time)

1. `lit_wiki_pending` writes `raw/compile_batch.json` for up to 25 papers and returns its path.
   `n_without_text` counts papers with no abstract; for those run `lit_curate` with
   `action=enrich` first, or leave them out.
2. Read the batch. For every item write one result using **only that item's** title, abstract
   and `fulltext_excerpt`:

   ```json
   {"results": [{
     "publication_id": "<copy exactly>",
     "summary": "2-4 sentences of synthesis, not the abstract's first lines",
     "key_claims": ["3-6 claims the text supports"],
     "methods": ["techniques named in the text"],
     "concepts": [{"name": "T cell exhaustion", "evidence": "short phrase from the text"}],
     "tags": ["short", "topic", "tags"],
     "confidence_notes": "say what the text does not settle"
   }]}
   ```

   Save it to the batch's `response_path`.
3. `lit_wiki_apply`. Each result is checked separately and the valid ones are written. Read
   `failed`, keep only those entries in the response file with their problems fixed, and call it
   again (a result that is applied twice just rewrites the same page).
4. Repeat from step 1 until `n_pending` is 0. Add `with_fulltext` when a PDF is saved and the
   abstract alone is thin; the page then records that it was written from full text.

### Grounding rules

- Use only the supplied text. No invented identifiers, sample sizes, p-values or results. If the
  text does not say it, leave it out or put it in `confidence_notes`.
- Describe findings; do not give clinical advice.
- Concept names: the thing the paper is about, in its usual form (a disease, a pathway, a method,
  a gene or protein symbol such as `LAYN`). Do not split a name at hyphens (`PD-L1` is one concept).
  Reuse an existing concept page's name when the paper is about the same thing.
- A paper marked `needs_review` gets a banner on its page. Do not blend its claims into concept
  definitions.

## Write concept definitions

A concept page starts with a placeholder Definition. When `lit_wiki_pending`, `lit_wiki_apply` or
`lit_wiki_lint` lists `concepts_needing_definition` (two or more papers mention it):

1. Open the concept page and the source pages it links.
2. Replace the placeholder under `## Definition` with 2-5 sentences using only what those pages
   say. Where papers disagree, say so and name them.
3. Optionally add an `## Across papers` section with one line per paper. Edit the file directly;
   no tool is needed.

## Organise the index

`wiki/index.md` is the table of contents. New papers land under `## Unfiled`.

1. Read the source pages' tags and concepts to propose 3-8 theme headings; show them to the user.
2. After they agree, move the bullets under the headings in `index.md` with the edit tool. A paper
   may appear under more than one heading. Keep every source page linked somewhere in the index.
3. `lit_wiki_lint` must then show no `missing_index_entry` errors.

## Check

`lit_wiki_lint` reports `missing_index_entry`, `broken_link`, `orphan_page`,
`stale_needs_review_banner`, `raw_without_source_page` and `concept_needs_definition`. Errors need
fixing; warnings and info are for the user's judgement.
