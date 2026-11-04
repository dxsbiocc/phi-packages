---
name: lit-enrich
description: >-
  Enrich ID-only publication shells with real title/abstract/authors/year/venue
  via PubMed efetch / Europe PMC (or --mock offline). Uses richest-merge; no
  invented fields. Prefer before wiki agent-compile when abstract is missing.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: enrich
      description: Fetch metadata for shell publications and merge into SQLite.
      run: [python, ./scripts/enrich.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          publication_id:
            type: string
          dry_run:
            type: boolean
          mock:
            type: boolean
          mock_fixture:
            type: string
            format: project-path
          force:
            type: boolean
          verify_oa:
            type: boolean
            description: Soft-verify candidate OA URLs after enrich; never store PDF bytes
          oa_mock_fixture:
            type: string
            format: project-path
        required:
          - project_root
          - collection_id
      approval: write
      output: ./schemas/enrich-result.json
      timeoutSeconds: 180
---

# lit-enrich (v0.5.1)

Fill metadata for **ID-only shells** (title empty/equals doi|pmid, abstract empty).

```bash
# After ids ingest (also runs automatically unless --no_enrich)
python scripts/enrich.py --project_root P --collection_id C

# Offline smoke
python scripts/enrich.py --project_root P --collection_id C --mock

# One pub / dry-run
python scripts/enrich.py --project_root P --collection_id C \
  --publication_id 'doi:10.1016/j.cell.2017.05.035' --dry-run
```

Fetch order: PMID → PubMed efetch; else DOI → PubMed `[doi]` then Europe PMC `DOI:"…"`.

工具名：`lit_enrich`。Wiki compile：若 abstract 仍空，先 enrich 再 `wiki_compile_request`。


## 0.5.1

`--verify_oa` stores `oa_status` / `oa_json` on the publication. Candidate URLs are never treated as owned PDFs.
