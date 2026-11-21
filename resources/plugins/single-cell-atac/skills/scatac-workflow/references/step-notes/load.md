Purpose: bring a scATAC input into a backed AnnData with per-cell fragment
metrics. Principle: snapatac2.import_fragments reads a 10x fragments.tsv.gz
(chrom, start, end, barcode, count) directly into a file-backed object and
computes n_fragment/frac_dup/frac_mito; a BAM is first converted with
make_fragment_file, a saved .h5ad is copied and re-backed, a 10x matrix is read
as peaks. Interpretation: the input kind is detected from the path, and the
genome's chrom_sizes define the tile space; read the reported n_fragment
distribution before filtering. Pitfalls: min_num_fragments drops barcodes below
the threshold AT IMPORT (an empty-droplet gate) — set it from the data, not a
tutorial default, or you can discard every cell; the working .h5ad is a temp
copy, never the user's input file; a fragments file needs the right chrom_sizes
or chromosomes mismatch. Raw-read alignment/Cell Ranger is upstream, not here.