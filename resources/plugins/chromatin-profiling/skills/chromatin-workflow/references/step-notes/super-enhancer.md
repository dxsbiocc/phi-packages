# Step note: ROSE super-enhancer calling (`super-enhancer`)

## Purpose
Stitch constituent enhancer peaks within `stitch_distance` (default 12.5 kb), integrate enhancer signal per domain, and separate **Super-Enhancers (SE)** from **Typical Enhancers (TE)** via the exact ROSE normalized tangent inflection point ($i^* = \arg\max_i (x_i - y_i)$).

## How it works
- Reads BED/narrowPeak/broadPeak constituent enhancer peaks (`signalValue` or `score` column).
- Optionally excludes promoter peaks overlapping TSS $\pm$ `tss_exclude_bp` (default 2,500 bp) when `tss_bed` or `gtf`/`txdb` is supplied.
- Stitches adjacent regions within `stitch_distance` using `GenomicRanges::reduce`, ranks summed signal ascending, and identifies the geometric slope = 1 inflection threshold.

## Reading the results
- `super_enhancers.tsv`: stitched domains ranked by integrated signal (`SE` vs `TE`, `n_constituents`, width).
- `super_enhancers.bed`: BED coordinates of called SEs only.
- `figures/super-enhancer-hockeystick.png`: ranked signal curve with cutoff lines.

## Pitfalls
- Plain 3-column BED files lack signal scores (the tool falls back to interval width weighting with a note); pass MACS2/SEACR narrowPeak/broadPeak or an explicit signal column for quantitative ROSE calling.
