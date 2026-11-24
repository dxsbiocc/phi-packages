# `mg_mag_traits` — AMR, virulence, and secondary-metabolite BGC summary

## Purpose

Aggregates specialized functional trait calls produced by upstream Nextflow
wrappers (`nf-core/amrfinderplus`, `nf-core/abricate`, `nf-core/rgi`, and
`nf-core/antismash`) into per-MAG summary tables (`mag-traits-summary.tsv`),
`bin × AMR class` and `bin × BGC product type` matrices, and an overview bar
chart.

## Inputs and attribution

- `--amrfinder`: one or more AMRFinderPlus TSV files (classifies `AMR`,
  `VIRULENCE`, and `STRESS` via `Element type`).
- `--abricate`: one or more ABRicate TSV files (uses `#FILE` column for bin
  attribution; `vfdb` / `virulence` databases count toward virulence hits).
- `--rgi`: one or more RGI/CARD TSV files (`Drug Class` counts).
- `--antismash`: per-bin antiSMASH `*.json` files or a directory of per-bin
  subdirectories containing antiSMASH JSON outputs (`region` features).
- `--min_identity` / `--min_coverage`: optional percent floors applied to
  AMRFinderPlus, ABRicate, and RGI hits (default `0`, preserving upstream calls).
- `--bin_map`: optional `file,bin` mapping table overriding filename attribution.

## Interpretation and pitfalls

- **Homology ≠ phenotype**: detected AMR or virulence genes indicate genetic
  potential, not confirmed clinical resistance or expression.
- **Draft MAG fragmentation**: antiSMASH BGC calls on fragmented contigs can be
  truncated at contig edges (`contig_edge` regions); interpret BGC counts
  alongside `mg_assembly_stats` N50 and `mg_bin_qc` completeness.
- **Deduplication across tools**: combining AMRFinderPlus, ABRicate, and RGI on
  the same MAGs will sum hits from each tool separately rather than harmonizing
  overlapping coordinates — report `sources` clearly or run `hamronization`
  upstream if strict cross-tool deduplication is required.
