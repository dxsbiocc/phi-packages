# Bin filtering

**Purpose.** Materialize the quality decision: a directory of passing bin
FASTAs plus a manifest, ready for dRep/GTDB-Tk/CoverM/annotation.

**Principle.** Filtering is a pure completeness/contamination gate on the
quality table; files are copied, never modified. Bin identity is the file
stem (extension stripped), matching CheckM2 `Name`, CheckM `Bin Id` and
GTDB-Tk `user_genome` conventions.

**Interpretation.** The pass rate is a pipeline health metric: healthy
complex-community runs often pass 30–70% of bins at MQ; near-zero means
upstream failure (depth, host contamination, wrong binner settings), not a
reason to relax thresholds.

**Pitfalls.** The default 50/10 floor is MIMAG MQ — appropriate for
community surveys, too lax for reference-genome claims (use 90/5). Copying
duplicates disk usage; on tight disks prefer symlinks via the project
environment instead. A quality table listing bins absent from bins_dir
fails the run — reconcile names before retrying, do not edit the table to
make it pass.