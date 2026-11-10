# Importing and exporting

## lit_import

| Input | What happens |
|--|--|
| `.bib`, `.ris` | Each entry becomes a paper. IDs in the entry are used for deduplication. |
| A text file of DOIs / PMIDs / PMCIDs / arXiv ids | One paper per id; the title and abstract are then looked up automatically. |
| `.jsonl` of hits | Same format `lit_search` writes. |
| `.pdf` | Text is read with `pdftotext`; a PDF with no text layer is rasterized with `pdftoppm` and read with `tesseract` (first 40 pages) when both are installed. The paper's own DOI, PMID, PMCID or arXiv id is taken from the first two pages (never the reference list) or the file name. The path and content hash are recorded; the PDF is not copied into SQLite. |
| A folder | Every `.pdf`, `.bib`, `.ris` and `.jsonl` inside, recursively. |

- With an id found, the title, abstract, year and authors are looked up from PubMed, Europe PMC
  or Crossref. The PDF's text is not stored as the abstract; read it with
  `lit_get with_fulltext`.
- A PDF with no recognisable id is saved under a temporary id (`tmp:...`) with a title guessed
  from the PDF metadata or first page. Enrichment then looks that title up in Crossref and
  accepts it only when the match is near-identical; the `tmp:` record is merged into the found
  paper (`matched_by_title` in the result: tell the user which papers were identified this way).
  If nothing matches, tell the user which files stay unidentified and ask for the DOI or PMID.
  Never guess one.
- `no_enrich` skips the metadata lookup (use it offline).
- `verify_oa` checks whether open-access links serve a PDF. Without a link in the record it asks
  OpenAlex for one by DOI. A link is never treated as an owned file.
- `access_origin` is a free label such as `institutional` or `author` for files the user
  obtained themselves. It is stored with the file; it does not fetch anything.
- Saved files stay where they are. If they move, use `lit_curate action=rebind`.

## lit_export

| `format` | Output |
|--|--|
| `bibtex` | `.bib`, entry types `article` or `misc`, double-braced titles so capitalisation survives |
| `ris` | `.ris` for EndNote, Zotero, Mendeley |
| `csl-json` | For Zotero, Pandoc and citeproc |
| `csv` | Spreadsheet view with the cite key, authors, venue, DOI, PMID, tags |
| `search_log` | `.jsonl` of every paper with the run that found it, for documenting a search |

Cite keys look like `smith2020immune` (first author, year, first title word of four or more letters) and get
a letter suffix on a clash. They are stable for the same set of papers but can change when
papers are added, so export again before citing rather than reusing an old key.

Papers with `needs_review` are left out unless `include_needs_review` is set; the result
reports how many were skipped.
