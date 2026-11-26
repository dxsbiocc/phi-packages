# Citation Faithfulness & BibTeX Audit Branch

Load this reference ONLY when auditing manuscript citations, DOIs/PMIDs, retractions, or `.bib` files.

## 1. Retrieval Routing & Offline Degradation

1. **`literature-library` plugin (`LiteratureLibrary` agent, if enabled)**: `lit_search` (sources: `pubmed`, `europepmc`, `openalex`, `biorxiv`) / `lit_get(..., with_fulltext=True)` (reads locally saved PDFs) for metadata and claim alignment; `lit_curate(action="versions", http=True)` for preprint -> published resolution.
2. **Database MCPs (if enabled)**: `pubmed`, `biorxiv` (or user-configured external literature MCPs).
3. **Local Workspace**: Cross-check `.tex`/`.md` keys against `.bib`/`.ris` and local PDFs.
4. **Offline Fallback**: If 1–2 are unavailable, run **Structural Audit Mode Only** (orphan `\cite{}`, unused `.bib`, duplicates, missing fields, malformed DOIs) and mark external veracity `[UNVERIFIED — RETRIEVAL TOOL REQUIRED]`.

## 2. Audit Dimensions & Severity

| Axis | What to Catch | Severity |
| --- | --- | --- |
| **Existence & Metadata** | Hallucinated titles/DOIs/PMIDs; wrong year, journal, or first author | `Critical` |
| **Retraction / Erratum** | Retracted articles or expressions of concern cited as valid evidence | `Critical` |
| **Claim-Citation Faithfulness** | Paper shows mouse correlation, manuscript cites it for human causality | `Major` |
| **Primary vs Secondary** | Citing a review article for a primary experimental number or cohort `n` | `Major` |
| **Preprint Supersession** | Citing an old preprint now published in a peer-reviewed journal | `Minor / Major` |
| **BibTeX / Key Hygiene** | Undefined `\cite{}`, unused `.bib` entries, duplicate keys, broken braces (`{RNA}`) | `Minor` |

## 3. Output Format

Deliver: `Audit Mode (Full vs Structural Only)` + Table (`Key | Status | Severity | Issue | Verified Metadata / Required Fix`) + cleaned BibTeX diff if requested.
