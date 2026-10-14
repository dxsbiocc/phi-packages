# Literature sources in this catalog

Connectors are MCP packages under `resources/connectors/`. They are not skills. Use whatever tools the installed connector exposes; names vary by server version, so read the tool list before the first call.

| Question | Connector id | What it covers | Auth in the package |
| --- | --- | --- | --- |
| Published articles and available PMC full text | `pubmed` | Biomedical literature, abstracts, PMC full text when the server returns it | none |
| Preprints | `biorxiv` | bioRxiv and medRxiv | none |
| Interventional studies | `clinical-trials` | ClinicalTrials.gov | none |

Related connectors, not literature search:

| Question | Connector id |
| --- | --- |
| Target–disease–drug associations | `open-targets` |
| Cancer studies, mutations, and clinical data in cBioPortal | `cbioportal` |

`pubmed`, `biorxiv`, `clinical-trials`, and `open-targets` declare `auth: none`. `cbioportal` declares `auth: oauth`. If a call fails for authentication, report that and do not replace the result with an unsourced citation.

## Search shape

- Put gene symbols, disease names, and method words in the query the connector accepts. Avoid a long natural-language sentence when the tool takes a keyword query.
- Apply the user's date bound, species, and article type when the tool has those fields.
- Read a small page of results, then narrow. Do not dump an unfiltered first page into the note.
- Prefer the latest indexed record when a preprint and a later PMID describe the same work. Keep both identifiers and say which one is the preprint.

## Identifiers

Copy these from the tool payload. Do not reconstruct them.

| Source | Keep |
| --- | --- |
| PubMed | PMID, and PMCID or DOI when present |
| Preprint | DOI and the server (bioRxiv or medRxiv) |
| Trial | NCT number |

A citation without one of those, taken from the same response, does not go in the note.
