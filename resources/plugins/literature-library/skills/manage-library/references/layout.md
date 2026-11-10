# Project layout

```
<project>/
  library.sqlite                 papers, aliases, collections, searches, files, version links
  exports/                       lit_export output
  collections/<collection_id>/
    runs/<run_id>/               run.json (query, sources, dates, counts) and hits.jsonl
    expands/<expand_id>/         hits and seed of a saved lit_expand
    wiki-vault/                  see skill compile-wiki
```

Main tables: `publications` (one row per paper), `aliases` (every id that points at a paper),
`collection_pubs`, `runs` and `run_hits` (what each search returned and what dedup did),
`pub_files` (saved file path, content hash, how it was obtained), `publication_versions`.

Papers are shared across collections: the same paper in two collections is one row.
