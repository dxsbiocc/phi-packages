---
name: lit-coverage
description: >-
  Coverage report for a run_id and/or collection: search hits vs unique ingested
  pubs, dedupe_action counts, needs_review, missing strong IDs, unenriched shells,
  gold-set recall, out-of-window leakage, declared-but-uncontacted sources.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: coverage
      description: JSON coverage report; optional markdown summary.
      run: [python, ./scripts/coverage.py]
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
          collection_id:
            type: string
          md:
            type: boolean
          md_out:
            type: string
            format: project-path
        required: []
      approval: read
      output: ./schemas/coverage-result.json
      timeoutSeconds: 60
---

# lit-coverage

```bash
python skills/lit-coverage/scripts/coverage.py \
  --project_root ./demo-project --run_id run-... --md

python skills/lit-coverage/scripts/coverage.py \
  --project_root ./demo-project --collection_id col-...
```

Requires `--run_id` and/or `--collection_id`. Tool name: `lit_coverage`.

For gate thresholds + `runs.status`, prefer `lit-run-close` / `lit-search --close_run`.

Also reports (0.5.2):

- **gold_set recall** (found/missed; `evidence_weak` if empty or &lt;5)
- **out_of_window_controls** leakage
- **declared-but-uncontacted sources** (always listed with reasons)

Coverage JSON/MD may include `hints[]` when gold_set recall < 1.0 or shells/uncontacted sources exist (0.5.3; honest, no fake recall).
