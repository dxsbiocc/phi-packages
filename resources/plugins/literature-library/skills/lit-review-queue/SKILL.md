---
name: lit-review-queue
description: >-
  List publications needing human identity adjudication only (soft fingerprint
  conflicts, conflicting strong IDs, empty unenriched shells, pending
  preprint↔journal version links). JSON stdout; optional markdown summary.
  Not full-text review.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: review_queue
      description: Diggable identity review queue (JSON; optional --md / --md_out).
      run: [python, ./scripts/review_queue.py]
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
          collection_id:
            type: string
          limit:
            type: integer
          no_shells:
            type: boolean
          no_versions:
            type: boolean
          md:
            type: boolean
          md_out:
            type: string
            format: project-path
          dig:
            type: boolean
          write_report:
            type: string
            format: project-path
        required: []
      approval: read
      output: ./schemas/review-queue-result.json
      timeoutSeconds: 60
---

# lit-review-queue

Identity-only human queue (not full-text review). Reason codes:

| code | meaning |
|--|--|
| `soft_fingerprint_conflict` | same title+first_author+year soft fingerprint, different strong IDs |
| `conflicting_strong_ids` | draft tied two existing pubs / `#conflict-` row |
| `empty_unenriched_shell` | ID shell with no real title/abstract |
| `pending_version_link` | `publication_versions.needs_confirm=1` |

```bash
python skills/lit-review-queue/scripts/review_queue.py \
  --project_root ./demo-project --collection_id col-... --md
```

```bash
# Schedulable dig — empty backlog exits 0 with silent minimal JSON
python skills/lit-review-queue/scripts/review_queue.py \
  --project_root P --collection_id C --dig --write_report /tmp/lit-dig.md
```

`--dig` empty → `silent:true`, `n:0` (exit 0). Non-empty items include
`suggested_next_action`. `--write_report` writes markdown **only when n>0**.
Do not create Grok Bot routines here — CLI is schedulable for the parent.
Tool name: `lit_review_queue`.
