Purpose: build the genome-wide tile matrix and select informative features for
the embedding, then remove doublets. Principle: a tile matrix bins the genome
uniformly (bin_size bp) so every region is represented regardless of whether a
peak was called; select_features keeps the most variable tiles and excludes the
blacklist; scrublet simulates doublets from the selected features and scores each
cell's doublet probability. Interpretation: tiles are for clustering (uniform
coverage), called peaks come later for differential accessibility; read the
doublet rate against the expected rate for the loaded cell number. Pitfalls:
bin_size is a biological resolution decision — large tiles for coarse domains,
small for fine structure; scrublet's expected_doublet_rate scales with cell
number (10x loads more cells → more doublets); a 10x peak-matrix input already
has features so the tile step is skipped. Doublet filtering changes cell counts,
so it is recorded as a QC gate.