# literature-library (Literature Library) · 0.7.0

Project-level literature library for Phi: search open sources, import your own files, deduplicate,
export citations, and compile saved papers into a wiki. No API keys, no embeddings.

A project-level literature library: search open sources, import your own files, deduplicate
automatically, export citations, and compile the papers already saved into a wiki.

## Components

| Component | Role |
|--|--|
| Agent `LiteratureLibrary` | Role and boundaries, about 40 lines; the workflows all live in the skills |
| Skill `manage-library` | Create the library, search, import, expand, query, export, maintain (9 tools) |
| Skill `compile-wiki` | Write papers into a wiki: summaries, key claims, concept pages, thematic index (using the 3 wiki tools in `manage-library`) |

`manage-library` ships all the scripts and `compile-wiki` is documentation only, so there is no
code dependency between the two skills.

## Tools (prefix `lit_`)

| Tool | Purpose |
|--|--|
| `lit_status` | Overview of the library and what to do next |
| `lit_collection` | Create or list collections (topics) |
| `lit_search` | Search **external** sources and save the hits with deduplication |
| `lit_import` | Import `.bib`, `.ris`, id lists, PDFs or a whole folder |
| `lit_expand` | From a seed paper, find similar papers, papers it cites and papers that cite it (preview only by default) |
| `lit_query` | Search the papers and wiki pages **inside** the library |
| `lit_get` | Read one paper's full record, files, versions and wiki page, optionally with PDF text |
| `lit_export` | Export BibTeX, RIS, CSL-JSON, CSV or the search log |
| `lit_curate` | Human review queue, metadata enrichment, preprint-to-published links, duplicate merging, rebinding files that moved |
| `lit_wiki_pending` | List papers that have no wiki page yet and write one compile batch |
| `lit_wiki_apply` | Validate the results the agent wrote and write them to disk |
| `lit_wiki_lint` | Health-check the wiki |

## Requirements

- Environment `phi:python@1`; the scripts use only the Python standard library, with no third-party
  dependencies. Network access is used only for `lit_search`, `lit_expand` and metadata enrichment.
- Reading text from a local PDF needs `pdftotext` (poppler); for a PDF without a text layer,
  `tesseract` or `LIT_PDF_OCR_CMD` are optional, and a missing one only produces a warning.
- Optional environment variables: `OPENALEX_MAILTO`, `CROSSREF_MAILTO`, `CONTACT_EMAIL` (polite
  contact addresses).

## Tests

```bash
python3 -m pytest resources/plugins/literature-library/tests -q
```

Fully offline. `tests/test_plugin_contract.py` checks the manifest, the skill frontmatter, and that
every tool argument is a real argument of the corresponding script.

## Boundaries

No citation fact-checking, evidence scoring, hypothesis work or omics analysis; no connection to
Google Scholar, CNKI or commercial literature gateways; no embedded MCP server, and no bundled keys.
Design trade-offs are in [DESIGN.md](DESIGN.md); version changes are in [CHANGELOG.md](CHANGELOG.md).