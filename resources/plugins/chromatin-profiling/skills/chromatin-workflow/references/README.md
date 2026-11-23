# references — chromatin-workflow

Organised by function, not piled up by source. Two kinds of reader, two entry points.

## Agent — how to run it, how to read it

- [`analysis-playbook.md`](analysis-playbook.md) — the operating manual. Each phase is a decision process:
  which tool to **call**, which JSON fields and figures to **read**, which parameters to **decide** from the
  numbers read, and which **red lines** to respect. The threshold starting points target ChIP/CUT&Tag data and
  are not defaults; read the distribution first, then adapt. The engines are R/Bioconductor: DiffBind 3.20
  (occupancy+affinity), DESeq2/edgeR (differential), csaw 1.44 (broad mark window method), ChIPseeker 1.46
  (annotation), ComplexHeatmap (UpSet), ggplot2 (PCA/MA). **Pure downstream**: starts from called
  peaks/count table/indexed BAM/bigWig and never touches FASTQ.
- [`red-line-cheatsheet.md`](red-line-cheatsheet.md) — hard red lines §1–§9 × common misreadings quick-reference
  table, for use at runtime; the full decision process is still governed by analysis-playbook.md.
- [`step-notes/`](step-notes/) — one short note per step (purpose/principle/interpretation/pitfalls), named by
  step id and paired with `assets/fragments/<step>.R`. These are exactly the teaching Markdown cells that the
  notebook renderer embeds.

## Maintainers — implementation and provenance

- [`maintainers/implementation-notes.md`](maintainers/implementation-notes.md) — what each step computes, how the
  fragment library/adapter/report layer divide the work, the DiffBind/csaw/ChIPseeker API mapping, the lean
  annotation and blacklist self-containment constraints, and the JSON contract.

## Sources and licences

The knowledge is an original distillation of the open sources below; no upstream text/figures/data are copied:

| Source | URL | Licence |
| --- | --- | --- |
| DiffBind vignette (Stark & Brown) | https://bioconductor.org/packages/DiffBind | Artistic-2.0 |
| csaw vignette (Lun & Smyth) | https://bioconductor.org/packages/csaw | Artistic-2.0 |
| ChIPseeker vignette (Yu et al.) | https://bioconductor.org/packages/ChIPseeker | Artistic-2.0 |
| ENCODE ChIP-seq guidelines (Landt et al.) | doi:10.1101/gr.136184.111 | Genome Res 22:1813 (2012) |
| ENCODE Blacklist (Amemiya, Kundaje, Boyle) | doi:10.1038/s41598-019-45839-z | Sci Rep 9:9354 (2019); Boyle-Lab/Blacklist ENCODE v2 |
| nf-core/chipseq (upstream wrapper) | https://github.com/nf-core/chipseq | MIT |

The ENCODE hg38/mm10 blacklist is inlined as base64-gzip in `assets/fragments/blacklist-filter.R`, so both the
script tools and the rendered notebook can filter artifact regions offline. The target-genome annotation
(GTF/GFF3 or TxDb) and org.db are supplied by the user via `input-path` or through a project environment, and are
never bundled for the target genome nor installed into the read-only plugin environment. The code is original to
this plugin and runs in the private `plugin:chip` environment.