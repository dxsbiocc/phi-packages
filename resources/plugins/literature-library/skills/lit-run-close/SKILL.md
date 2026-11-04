---
name: lit-run-close
description: >-
  Close a frozen search run: write coverage.json (+ md), evaluate coverage gate
  thresholds, mark runs.status passed|needs_attention. Schedulable CLI; no
  commercial sources.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: run_close
      description: Coverage gate + mark run status after search/ingest.
      run: [python, ./scripts/run_close.py]
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
          run_id:
            type: string
          min_unique_pubs:
            type: integer
          max_needs_review:
            type: integer
          max_version_pending:
            type: integer
          max_unenriched_shells:
            type: integer
          fail_under_coverage:
            type: boolean
          coverage_json:
            type: string
            format: project-path
          coverage_md:
            type: string
            format: project-path
        required:
          - run_id
      approval: write
      output: ./schemas/run-close-result.json
      timeoutSeconds: 60
---

# lit-run-close

After `lit-run-config` → `lit-search` (+ ingest), close the run:

```bash
python skills/lit-run-close/scripts/run_close.py \
  --project_root P --run_id run-... \
  --min_unique_pubs 5 --max_needs_review 20

# Exit 2 when gate fails:
python skills/lit-run-close/scripts/run_close.py \
  --project_root P --run_id run-... \
  --min_unique_pubs 10 --fail_under_coverage
```

Writes `runs/<id>/coverage.json` + `coverage.md`, updates `run.json` / `runs.status`
to `passed` or `needs_attention`. Same thresholds available on `lit-search --close_run`.

Tool name: `lit_run_close`.
