Purpose: coverage profile around the TSS and a per-peak heatmap — the shape of the
signal, not just its location. Principle: plotAvgProf2 takes the peaks and the user
TxDb (upstream/downstream), builds a tag matrix internally (getTagMatrix) and
averages it with a 95% confidence band; peakHeatmap draws per-peak rows. A
sharp TF peak concentrates at the motif; a broad mark spreads over the gene body.
Interpretation: the average TSS profile shows where signal piles up relative to
transcription start; the heatmap shows per-peak heterogeneity. bigWig tracks, when
present, give the same view genome-wide. Pitfalls: needs a user TxDb (gtf_path or
txdb_path) — run annotate first or pass the annotation source; without one it
stops. It uses the consensus peaks when available, else the per-sample peak list.
The profile is descriptive signal shape, not a quantitative binding measure.