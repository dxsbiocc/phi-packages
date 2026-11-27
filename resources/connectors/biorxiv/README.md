bioRxiv & medRxiv provides unified access to Cold Spring Harbor Laboratory's life-science (`bioRxiv`) and clinical/health-science (`medRxiv`) preprint servers, combining indexed keyword search via Europe PMC with the official `api.biorxiv.org` date-interval, version-history, and journal-publication APIs.

### Capabilities

- Search bioRxiv and medRxiv preprints by keyword, gene, disease, method, or author with optional server, category, and date filters.
- Browse newly posted or revised preprints over any date window and subject category (`bioinformatics`, `genomics`, `cancer biology`, `immunology`, `neuroscience`, and others).
- Inspect full preprint metadata, abstracts, revision history (`v1` through latest version), corresponding institution, license, and direct PDF / JATS XML links.
- Track whether a preprint has subsequently been published in a peer-reviewed journal and resolve its journal DOI (`published_as`).

### Research & Workflow Use Cases

- Monitor emerging computational methods, single-cell atlases, and clinical preprints ahead of journal publication.
- Check whether a cited `10.1101/...` preprint now has a peer-reviewed journal version to cite in manuscripts.
- Compare revision dates and version counts when evaluating rapidly evolving preprints.

### Tools

| Tool | Description |
| --- | --- |
| `usage_terms` | Return usage terms, preprint licensing notes, and routing details for bioRxiv and medRxiv. |
| `biorxiv_search_preprints` | Search bioRxiv and medRxiv preprints by topic, gene, disease, or author with server, category, and date filters. |
| `biorxiv_browse_recent` | Browse newly posted or revised preprints on bioRxiv or medRxiv over a date window via CSHL's official details API. |
| `biorxiv_get_preprint` | Retrieve full preprint metadata, abstract, version history, PDF/XML links, and linked peer-reviewed journal DOI. |
| `biorxiv_get_published` | Check peer-reviewed journal publication status for a preprint DOI or browse recent preprint-to-journal publications. |

### Usage Notes

- No sign-in or API key is required. Hosted read-only by [Phi Packages](https://github.com/dxsbiocc/phi-packages).
- Preprints have not undergone peer review; always verify clinical or experimental claims and check `published_as` for peer-reviewed updates.
- Respect each preprint's Creative Commons or author-selected license on [bioRxiv](https://www.biorxiv.org/) and [medRxiv](https://www.medrxiv.org/).
