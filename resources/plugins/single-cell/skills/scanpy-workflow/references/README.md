# References index

Grouped by function: agent decision docs first, maintainer provenance under
`maintainers/`.

## Agent — how to run and interpret

- [analysis-playbook.md](analysis-playbook.md) — **read this when analyzing**:
  per-stage decision procedures (call → read → decide → red lines), threshold
  starting points, and the JSON field guide.
- [step-notes/](step-notes/) — one note per analysis step (purpose / principle /
  interpretation / pitfalls), paired to the fragment by step id. These are the
  notebook's teaching cells; read the note beside the fragment you run.
- [cross-route-handoff.md](cross-route-handoff.md) — moving objects between the
  R and Python routes (pseudobulk tables or explicit conversion).

## Maintainers — implementation and provenance

- [maintainers/implementation-notes.md](maintainers/implementation-notes.md) —
  what each script tool computes and how the NBIS labs map onto the chain;
  read when adapting a step.
- [maintainers/nbis-porting.md](maintainers/nbis-porting.md) — porting
  decisions, provenance and license (CC-BY-4.0 per the repository LICENSE).