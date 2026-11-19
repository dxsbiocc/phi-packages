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

## Maintainers — implementation and provenance

- [maintainers/implementation-notes.md](maintainers/implementation-notes.md) —
  what each script tool computes and how the scNotebooks modules map onto the
  chain; read when adapting a step.
- [maintainers/scnotebooks-curriculum.md](maintainers/scnotebooks-curriculum.md) —
  the 13-module teaching curriculum index.
- [maintainers/scnotebooks-module-map.md](maintainers/scnotebooks-module-map.md) —
  per-module package facts and runnability.
- [maintainers/scnotebooks-porting.md](maintainers/scnotebooks-porting.md) —
  porting decisions and provenance.