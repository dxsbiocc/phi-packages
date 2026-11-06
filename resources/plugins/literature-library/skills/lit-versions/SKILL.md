---
name: lit-versions
description: >-
  Preprint↔journal version chain: detect candidates (offline heuristics +
  optional OpenAlex/Crossref HTTP), propose publication_versions links, confirm
  (+ optional merge into journal DOI canonical; preprint IDs kept as aliases).
  Wiki sources repoint on merge. No commercial sources.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: versions
      description: detect | apply | propose | confirm version links.
      run: [python, ./scripts/versions.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          db_path:
            type: string
            format: project-path
          action:
            type: string
            description: detect|apply|propose|confirm
          collection_id:
            type: string
          from_id:
            type: string
          to_id:
            type: string
          relation:
            type: string
          confidence:
            type: number
          link_id:
            type: integer
          no_merge:
            type: boolean
          auto_merge_high_confidence:
            type: boolean
          evidence_json:
            type: string
          http:
            type: boolean
          providers:
            type: string
          mock:
            type: boolean
          fixture_dir:
            type: string
          no_wiki_repoint:
            type: boolean
          preview_fetch:
            type: boolean
        required: []
      approval: write
      output: ./schemas/versions-result.json
      timeoutSeconds: 180
---

# lit-versions

Schema table `publication_versions`: `from_id`, `to_id`, `relation`
(`preprint_of` / `published_as` / `same_work`), `confidence`, `needs_confirm`, `evidence_json`.

```bash
# dry detect (offline heuristics)
python skills/lit-versions/scripts/versions.py --project_root P --action detect

# offline HTTP fixtures (OpenAlex locations + Crossref relation)
python skills/lit-versions/scripts/versions.py --project_root P --action detect --mock

# live OpenAlex / Crossref (soft-fail; public HTTP only)
python skills/lit-versions/scripts/versions.py --project_root P --action apply \
  --http --auto_merge_high_confidence

# write pending links
python skills/lit-versions/scripts/versions.py --project_root P --action apply

# high-confidence → merge into journal DOI (preprint DOI/arxiv remain aliases);
# wiki/sources/<old_safe>.md repointed (+ stub) unless --no_wiki_repoint
python skills/lit-versions/scripts/versions.py --project_root P --action apply \
  --auto_merge_high_confidence

# manual propose / confirm
python skills/lit-versions/scripts/versions.py --project_root P --action propose \
  --from_id arxiv:1234.5678 --to_id doi:10.1038/s41586-xxx --confidence 0.8
python skills/lit-versions/scripts/versions.py --project_root P --action confirm --link_id 1
```

Detection: **offline** title/soft-fingerprint/year, plus optional **HTTP**
(OpenAlex multi-location DOI; Crossref `relation`).

```bash
# Missing opposite DOI → metadata-only shell + pending link (no PDF)
python skills/lit-versions/scripts/versions.py --project_root P \
  --action apply --mock --preview_fetch
python skills/lit-versions/scripts/versions.py --project_root P \
  --action confirm --link_id 1   # merge to journal DOI preferred
```

No SmartLib / commercial gateways. Tool name: `lit_versions`.
