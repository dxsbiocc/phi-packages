# Vault layout and page formats

```
wiki-vault/
  AGENTS.md       conventions (written once; do not delete)
  vault.json      schema label
  raw/
    publications/<safe_id>/{meta.json, body.txt}   what each page was written from
    compile_batch.json, compile_batch_response.json   the current batch
  wiki/
    index.md      table of contents (themes)
    log.md        append-only log; one line per write
    sources/<safe_id>.md
    concepts/<slug>.md
```

`safe_id` is the paper's id with `:` and `/` replaced by `_`, for example
`doi_10.1016_j.cell.2017.05.035`. `slug` is the lower-cased concept name with other characters
turned into `-`.

## Source page

Header bullets (publication id, year, venue, `needs_review`, what it was written from, raw
snapshot), then `## Identifiers`, `## Summary`, `## Key claims`, `## Methods`, `## Abstract`,
`## Concepts` (links to concept pages), optional `## Confidence notes`, and `## Related sources`
where you can link other pages by hand.

## Concept page

```
# <Name>

Slug: `<slug>`

## Definition
<placeholder until you write it>

## Sources
- [<paper title>](../sources/<safe_id>.md) — <evidence phrase>

## Related concepts
```

Paper bullets under `## Sources` are added automatically. Everything else is yours to edit;
keep the `## Sources` heading so later papers can be appended.

## What you may edit

`wiki/` pages are editable markdown: prefer appending and linking over rewriting. `raw/` is a
record of what pages were written from: add, never rewrite. Merging two paper ids renames the
source page and leaves a redirect stub at the old name.
