# references — scatac-workflow

Organized by function, not piled up by source. Two kinds of reader, two entry points.

## Agent — how to run it, how to interpret it

- [`analysis-playbook.md`](analysis-playbook.md) — the operating manual. Each stage is one decision process:
  which tool to **call**, which JSON fields and figures to **read**, what to **decide** from the numbers read, and which **red lines** to respect.
  The threshold starting points are tuned for scATAC data, not defaults; read the distributions first, then adapt. The engine is snapatac2 2.10
  (file-backed, Rust); scanpy 1.12 is used only for the gene-activity downstream and matplotlib panels.
- [`red-line-cheatsheet.md`](red-line-cheatsheet.md) — §1–§10 hard red lines × common misreadings quick-reference table,
  for use at runtime; the full decision process still follows analysis-playbook.md.
- [`step-notes/`](step-notes/) — one short note per step (purpose/principle/interpretation/pitfalls), named by step id,
  paired with `assets/fragments/<step>.py`. These are exactly the teaching Markdown cells the notebook renderer embeds.

## Maintainers — implementation and provenance

- [`maintainers/implementation-notes.md`](maintainers/implementation-notes.md) — what each step computes,
  how the fragment library/adapters/report layer divide the work, the snapatac2 API mapping, the backed-object constraints and the JSON contract.

## Sources and licenses

The knowledge is an original distillation of the following open sources; no upstream text/figures/data are copied:

| Source | URL | License |
| --- | --- | --- |
| SnapATAC2 documentation and tutorials (Zhang et al.) | https://github.com/kaizhang/SnapATAC2 | MIT |
| Single-cell best practices (chromatin accessibility chapter) | https://www.sc-best-practices.org/ | Apache 2.0; Heumos et al., Nat Rev Genet 2023 |
| ENCODE Blacklist (Amemiya, Kundaje, Boyle) | doi:10.1038/s41598-019-45839-z | Sci Rep 9, 9354 (2019); Boyle-Lab/Blacklist ENCODE v2 |

The ENCODE hg38/mm10 blacklist is inlined as base64-gzip in `assets/fragments/blacklist.py`,
so both the script tools and the rendered notebook can filter artifact regions offline. snapatac2 bundles
chrom_sizes for hg38/mm10/hg19/mm39 offline; the gene annotation (GFF3) and genome FASTA are lazy downloads,
supplied by the user through `input-path` or downloaded once into the project directory, and are never installed into the read-only plugin environment.
The code is original to this plugin and runs in the private `plugin:scatac` environment.