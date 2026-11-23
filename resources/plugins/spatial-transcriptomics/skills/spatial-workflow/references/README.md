# references — spatial-workflow

Organised by function, not piled up by source. Two kinds of reader, two entry points.

## Agent — how to run it, how to read it

- [`analysis-playbook.md`](analysis-playbook.md) — the operating manual. Each stage is a decision
  procedure: **call** which tool, **read** which JSON fields and figures, **decide** the parameters
  from the numbers you read, and respect the **red lines**. The threshold starting points target
  Visium/Stereo-seq data and are not defaults; read the distribution first, then adapt.
- [`red-line-cheatsheet.md`](red-line-cheatsheet.md) — a quick-reference table of the §1–§10 hard
  red lines × common misreadings, for checking against at runtime; the full decision procedure is
  still analysis-playbook.md.
- [`step-notes/`](step-notes/) — one short note per step (purpose/principle/interpretation/pitfalls),
  named by step id and paired with `assets/fragments/<step>.py`. These are exactly the teaching
  Markdown cells the notebook renderer embeds.

## Maintainers — implementation and provenance

- [`maintainers/implementation-notes.md`](maintainers/implementation-notes.md) — what each step
  computes, how the fragment library/adapter/report layer divide the work, and the JSON contract.
- [`maintainers/osta-map.md`](maintainers/osta-map.md) — the OSTA (Bioconductor) chapter to
  plugin-Python mapping, and how R concepts translate to scanpy/squidpy.
- [`maintainers/stereo-seq-methods.md`](maintainers/stereo-seq-methods.md) — Stereo-seq data formats
  (GEM/GEF, square bin/cellbin) and the methodological points of the 7 application papers.

## Sources and licences

The knowledge is an original distillation of the following open textbooks; no upstream
text/figures/data are copied:

| Source | URL | Licence |
| --- | --- | --- |
| OSTA (Orchestrating Spatial Transcriptomics Analysis with Bioconductor) | https://bioconductor.org/books/OSTA/ | CC BY 4.0 (pending final check) |
| Single-cell best practices (spatial chapters) | https://www.sc-best-practices.org/ | Apache 2.0; Heumos et al., Nat Rev Genet 2023 |

The Stereo-seq methodology follows the 7 papers provided by the user (Dev Cell 2022, Zool Res 2024,
Cell Discov 2024, Nat Genet 2024, Cell Rep 2024×2, Cell Metab 2024); see
maintainers/stereo-seq-methods.md.