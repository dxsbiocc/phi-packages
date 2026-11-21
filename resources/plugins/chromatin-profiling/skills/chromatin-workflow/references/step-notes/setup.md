Purpose: load the R/Bioconductor stack and initialize the bookkeeping every later
step writes into. Principle: a single `setup` fragment is the common preamble for
both the script adapters and the rendered notebook, so the two paths share one
analysis source; `need()` fails loudly with a "use the plugin:chip environment"
hint when a core package is absent, while csaw is optional (a platform gap, not a
hard requirement). Interpretation: it creates the shared variables — peak_list,
consensus_peaks, occupancy_matrix, count_matrix, normalized_counts, dba_object,
csaw_result, diff_results, qc_metrics, blacklist_removed, annotation_result, txdb,
sample_sheet, pca_result, figure_paths, analysis_steps — and defines show_plot,
which prints inline only under IRkernel and otherwise writes a ragg PNG into
figures_dir. Pitfalls: figures are PNG (ggplot2/ChIPseeker/ComplexHeatmap), not
plotly; a missing core package aborts here rather than mid-analysis, so run
chip_env_probe first to know what the managed env actually ships.