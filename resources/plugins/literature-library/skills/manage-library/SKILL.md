---
name: manage-library
description: >-
  Build and maintain a project-level literature library: search PubMed, Europe
  PMC, OpenAlex, Crossref, arXiv and bioRxiv/medRxiv, import BibTeX, RIS, ID
  lists and PDF files, deduplicate, query saved papers, follow citations from a
  seed paper, and export BibTeX, RIS or CSL-JSON. Use for 找文献, 检索文献,
  文献库, 导入 bib/ris/PDF, 去重, 文献整理, 导出参考文献, 引文扩展, 综述素材.
  Not for omics analysis, citation fact-checking, or one-off lookups of a single
  paper (the PubMed connector does that).
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: status
      description: >-
        Overview of the project's library and what to do next: paper counts,
        papers needing review or an abstract, missing files, wiki pages still to
        write, recent searches. Call this first in a new session or after a long
        gap. Reads only.
      run: [python, ./scripts/status.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: input-path
          collection_id:
            type: string
            description: Limit to one collection; default all
        required: [project_root]
      approval: read
      timeoutSeconds: 60
    - name: collection
      description: >-
        Create a collection (a topic container with its own wiki vault) or list
        collections with paper counts. The project gets its library.sqlite on
        the first create.
      run: [python, ./scripts/collection.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          action:
            type: string
            description: create or list
          name:
            type: string
            description: Topic name, required for create
          description:
            type: string
          no_wiki:
            type: boolean
            description: Do not create the wiki vault
        required: [project_root, action]
      approval: write
      timeoutSeconds: 60
    - name: search
      description: >-
        Search EXTERNAL literature sources and save the hits into the library
        with deduplication. One call records the search, fetches hits per
        source and ingests them. To search what is already saved use lit_query.
      run: [python, ./scripts/search.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
            description: Omit when the project has exactly one collection
          query:
            type: string
            description: Search words; PubMed field tags such as [ti] and MeSH work for pubmed
          sources:
            type: string
            description: "Comma-separated, at most 4 of pubmed, europepmc, openalex, crossref, arxiv, biorxiv, medrxiv. Default pubmed,europepmc"
          date_from:
            type: string
            description: YYYY or YYYY-MM-DD
          date_to:
            type: string
          limit:
            type: integer
            description: Hits per source, 1-50, default 20
          skip_ingest:
            type: boolean
            description: Fetch and write hits.jsonl only; do not save into the library
        required: [project_root, query]
      approval: write
      timeoutSeconds: 180
    - name: import
      description: >-
        Save papers from files you already have: .bib, .ris, a DOI/PMID list,
        a .jsonl of hits, a PDF, or a folder of them. Deduplicates, records the
        file path and content hash, and fills missing titles and abstracts.
      run: [python, ./scripts/ingest.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          input_path:
            type: string
            format: input-path
            description: A file or a folder
          format:
            type: string
            description: bib, ris, ids, pdf, pdf_text or hits; detected when empty
          access_origin:
            type: string
            description: How the file was obtained when it is not a plain download, e.g. institutional. Recorded only.
          no_enrich:
            type: boolean
            description: Do not look up missing metadata afterwards
          verify_oa:
            type: boolean
            description: Check whether candidate open-access links really serve a PDF
        required: [project_root, input_path]
      approval: write
      timeoutSeconds: 300
    - name: expand
      description: >-
        Find papers related to one saved paper: similar, citing it, or cited by
        it. One round per call, and the hits are only previewed unless
        confirm_write is set. Follow citations further by calling again with a
        hit as the new seed.
      run: [python, ./scripts/expand.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          seed_id:
            type: string
            description: Library id, DOI, PMID, PMCID, arXiv id or OpenAlex W id
          relation:
            type: string
            description: similar, citations or references
          limit:
            type: integer
            description: 1-50, default 10
          confirm_write:
            type: boolean
            description: Save the hits into the collection (default is preview only)
        required: [project_root, seed_id]
      approval: write
      timeoutSeconds: 180
    - name: query
      description: >-
        Search papers ALREADY SAVED in the library (title, abstract, tags) and
        their wiki pages. Several words must all match. Never calls an
        external source; use lit_search for that.
      run: [python, ./scripts/query.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: input-path
          collection_id:
            type: string
          query:
            type: string
          tag:
            type: string
          needs_review:
            type: string
            description: "0 or 1; omit for both"
          limit:
            type: integer
          no_wiki:
            type: boolean
        required: [project_root]
      approval: read
      timeoutSeconds: 60
    - name: get
      description: >-
        Everything saved about one paper: full record and abstract, tags,
        review reasons, saved files, preprint/journal links and wiki page.
        Optionally includes text read from its saved PDF.
      run: [python, ./scripts/get.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: input-path
          publication_id:
            type: string
            description: Library id, DOI, PMID, PMCID, arXiv id or OpenAlex id
          with_fulltext:
            type: boolean
          fulltext_chars:
            type: integer
          with_wiki:
            type: boolean
        required: [project_root, publication_id]
      approval: read
      timeoutSeconds: 120
    - name: export
      description: >-
        Write the collection as bibtex, ris, csl-json or csv (with stable cite
        keys), or as a search_log.jsonl record of what was found. Papers with
        unresolved identity are left out unless include_needs_review is set.
      run: [python, ./scripts/export.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          format:
            type: string
            description: bibtex (default), ris, csl-json, csv or search_log
          out_path:
            type: string
            format: project-path
            description: Default is exports/<collection>.<ext> in the project
          publication_ids:
            type: array
            items:
              type: string
            description: Export only these papers
          include_needs_review:
            type: boolean
        required: [project_root]
      approval: write
      timeoutSeconds: 60
    - name: curate
      description: >-
        Maintain what is saved. action=queue lists papers needing a human
        decision; enrich fills title, abstract, year and authors from the DOI, PMID, arXiv id or, failing those, an exact title match;
        versions finds or confirms preprint-journal links; merge_ids merges
        duplicates that share complementary IDs (dry run unless apply);
        rebind re-attaches files that moved, matched by content hash.
      run: [python, ./scripts/curate.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          action:
            type: string
            description: queue, enrich, versions, merge_ids or rebind
          collection_id:
            type: string
          publication_id:
            type: string
            description: enrich only this paper
          limit:
            type: integer
          versions_action:
            type: string
            description: detect (default), apply, propose or confirm
          from_id:
            type: string
          to_id:
            type: string
          link_id:
            type: integer
          http:
            type: boolean
            description: versions - also ask OpenAlex and Crossref
          auto_merge_high_confidence:
            type: boolean
          apply:
            type: boolean
            description: merge_ids - write the merges
          scan_dir:
            type: string
            format: input-path
            description: rebind - folder that now holds the files
          confirm:
            type: boolean
            description: rebind - write the new paths
          force:
            type: boolean
            description: enrich - look papers up again even if already looked up
        required: [project_root, action]
      approval: write
      timeoutSeconds: 300
    - name: wiki_pending
      description: >-
        List saved papers that have no wiki page yet and write one compile
        batch for the agent to fill (see skill compile-wiki). Also lists
        concept pages that still need a Definition.
      run: [python, ./scripts/wiki_pending.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          limit:
            type: integer
            description: Papers per batch, 1-25, default 10
          publication_ids:
            type: array
            items:
              type: string
            description: Rewrite these papers even if a page exists
          with_fulltext:
            type: boolean
            description: Add text from a saved PDF when there is one
          fulltext_chars:
            type: integer
        required: [project_root]
      approval: write
      timeoutSeconds: 120
    - name: wiki_apply
      description: >-
        Validate the agent's batch results and write source pages, concept
        pages, index links and tags. Each result is checked on its own.
      run: [python, ./scripts/wiki_apply.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          response:
            type: string
            format: input-path
            description: Default is raw/compile_batch_response.json in the vault
        required: [project_root]
      approval: write
      timeoutSeconds: 120
    - name: wiki_lint
      description: >-
        Check the wiki vault: papers missing from the index, broken links,
        orphan pages, stale review banners, concepts without a Definition.
      run: [python, ./scripts/wiki_lint.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: input-path
          collection_id:
            type: string
        required: [project_root]
      approval: read
      timeoutSeconds: 60
---

# Manage Library

A project keeps one `library.sqlite` and any number of **collections** (topics). Each
collection has its own search history and a wiki vault. SQLite is the source of truth for
which papers exist and who they are; the wiki is the readable layer on top (see
`skill://compile-wiki`).

Tools are `lit_*`. Search outside the library with `lit_search`; search inside it with
`lit_query`. They are not interchangeable.

## Start

Call `lit_status` when you begin or return to a project. Its `next_steps` say what is
waiting. If the project has no library, create a collection first (`lit_collection`
with `action=create`). With a single collection you can leave `collection_id` out
everywhere.

## Workflows

**Find papers on a topic.** `lit_search` with the query and 1-3 sources. Biomedical topics:
`pubmed,europepmc`; add `openalex` for breadth, `biorxiv` or `arxiv` for preprints. Read
`per_source` in the result: a source with `ok: false` failed and the others still ran. Report
inserted / merged / needs_review counts. Details and query tips: `references/search.md`.

**Save papers you already have.** `lit_import` with a file or folder. A DOI/PMID list works
and gets its metadata filled automatically. PDFs are matched to papers by the DOI or PMID in
their text or name. Formats and rules: `references/import-export.md`.

**Grow from a key paper.** `lit_expand` with `seed_id` and `relation` (`similar`, `citations`,
`references`). Show the preview to the user, then repeat with `confirm_write` for the ones
worth keeping. For a second round, use one of the new hits as the next seed.

**Read what is saved.** `lit_query` to find, `lit_get` to read one paper in full (add
`with_fulltext` when a PDF is saved and the abstract is not enough).

**Cite or hand off.** `lit_export` with `bibtex`, `ris`, `csl-json` or `csv`.

**Clean up.** `lit_curate`: `queue` shows what needs a human decision, `enrich` fills ID-only
papers, `versions` links preprints to their published versions, `merge_ids` and `rebind` repair
identity and moved files. Rules: `references/identity.md`.

**Write the wiki.** Switch to `skill://compile-wiki`.

## Rules

- Never invent identifiers, titles or abstracts. If something is missing, say so and use
  `lit_curate` with `action=enrich`.
- A paper marked `needs_review` has an uncertain identity. Do not merge, cite or summarise it
  as settled until the user decides.
- `lit_search` reaches the network. Keep to the sources the user needs and the default limit
  unless asked for more.
- Only import files the user is entitled to hold. `access_origin` records where a file came
  from; nothing here fetches paywalled text or handles credentials.
- Report absolute paths of files you created and the counts the tools returned.

## Boundaries

This library finds, stores, deduplicates, searches, exports and writes up papers. It does not
verify citations, score evidence, or analyse omics data. Other sources such as Google Scholar
or commercial literature gateways are not used.

Files live under the project: `library.sqlite`, `collections/<id>/runs/`,
`collections/<id>/wiki-vault/`, `exports/`. See `references/layout.md`.
