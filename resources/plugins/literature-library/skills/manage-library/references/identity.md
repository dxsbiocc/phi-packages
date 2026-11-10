# Identity and cleanup

SQLite decides which records are the same paper.

## How papers are matched

A paper's id is its strongest identifier: DOI, then PMID, PMCID, arXiv id, OpenAlex id.

- A new record that shares a strong id with a saved paper is merged into it and fills its
  gaps.
- A match on title, first author and year alone is not enough. Both records are kept and
  marked `needs_review` with the reason `soft_fingerprint_conflict`.
- Two saved papers whose strong ids contradict each other are marked
  `conflicting_strong_ids`.
- A record that has only a DOI and another that has only the PMID of the same article can be
  merged when the evidence is strong; `lit_curate action=merge_ids` shows the plan (dry run)
  and writes it only with `apply`.

## lit_curate actions

| Action | Use |
|--|--|
| `queue` | The list of papers a human has to decide about, with a reason code each: `soft_fingerprint_conflict`, `conflicting_strong_ids`, `empty_unenriched_shell`, `pending_version_link`. Ask the user; do not resolve soft conflicts yourself. |
| `enrich` | Fill title, abstract, year and authors from the DOI, PMID or arXiv id (PubMed, Europe PMC, Crossref, arXiv). A paper with no id but a full title is matched only when Crossref returns a near-identical title (`matched_by_title` in the result; tell the user) and a `tmp:` record then takes the found identity. A paper already looked up is skipped unless `force`. Needs the network. |
| `versions` | Preprint and journal versions of the same work. `detect` lists candidates, `apply` saves them as pending links, `confirm` with a `link_id` accepts one and merges the preprint into the journal record (the preprint id stays as an alias). `http` also asks OpenAlex and Crossref. `propose` records a pair you know about. |
| `merge_ids` | Merge records that are the same paper by complementary ids. Dry run unless `apply`. |
| `rebind` | After moving files, point the saved paths at the new locations by content hash. Dry run unless `confirm`. |

Merging changes a paper's id. The wiki pages follow automatically (the page is renamed and a
redirect stub is left at the old name).
