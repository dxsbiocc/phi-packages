Purpose: build a consensus peakset — the union of reproducible binding sites that
defines the region universe every sample is measured over. Principle: disjoin all
samples' peaks into non-overlapping intervals (GenomicRanges::disjoin), keep
intervals present in at least min_occupancy samples, then merge adjacent kept
intervals (reduce). Interpretation: n_consensus and occupancy_per_sample show how
many sites survived and how many each sample carries; the effective occupancy is
capped at the number of samples. A high min_occupancy with few replicates throws
away real but variable sites; a low one with many samples admits noise — choose
from replicate count and biology, and say why. Pitfalls: blacklist filtering must
run first or artifact regions enter the consensus. The consensus coordinates are
what the count matrix and differential test use, so this decision shapes every
downstream result; it is not a default to accept blindly.