#!/usr/bin/env Rscript

# Template-ID: bar-enrichment-dot
#
# Purpose:
#   Draw a ranked enrichment dotplot: terms on y, gene ratio on x,
#   point size = gene count, fill = -log10(adjusted p).
#
# Inputs:
#   One row per enriched term. Default example:
#     - Description: term name
#     - GeneRatio: gene ratio (already a fraction)
#     - p.adjust: adjusted p-value
#     - Count: gene count
#     - ONTOLOGY: optional facet / colour group (BP, CC, MF, KEGG)
#
# Output:
#   A PDF, PNG, or SVG enrichment dotplot.
#
# Dependencies:
#   ggplot2, ggprism, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, sort, facet, palette).
#   Extra ontologies or databases are extra rows, not extra ids.
#   A Change / up-down column is an optional colour group, not a new
#   id. Do not run enricher, GSEA, or clusterProfiler in this script.
#
# Scientific assumptions:
#   GeneRatio, Count, and p.adjust are supplied enrichment results.
#   -log10(p.adjust) is a display transform, not a new test. This
#   script does not compute ORA, GSEA, or gene-set overlap.

local({
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (!length(file_arg)) {
        stop("Run this file with Rscript.", call. = FALSE)
    }
    dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
    for (i in seq_len(8)) {
        candidate <- file.path(dir, "lib", "common.R")
        if (file.exists(candidate)) {
            source(normalizePath(candidate), chdir = FALSE)
            return(invisible())
        }
        parent <- dirname(dir)
        if (identical(parent, dir)) break
        dir <- parent
    }
    stop("Cannot find scripts/lib/common.R", call. = FALSE)
})

io <- parse_io_args()

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        description = "Description",
        gene_ratio = "GeneRatio",
        pvalue = "p.adjust",
        count = "Count",
        group = "ONTOLOGY"
    ),
    sort = "gene_ratio",
    facet = FALSE,
    palettes = list(
        fill = "Quantitative.BluGrn"
    ),
    size_range = c(2.8, 8.5),
    labels = list(
        x = "Gene ratio",
        fill = "-log10(p.adjust)",
        size = "Count"
    ),
    size = list(
        width = 7.6,
        height = 6.4
    )
)

load_packages(c("ggplot2", "ggprism", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
need <- c(
    config$columns$description,
    config$columns$gene_ratio,
    config$columns$pvalue,
    config$columns$count
)
require_columns(df, need)
df <- as.data.frame(df, stringsAsFactors = FALSE)

desc_col <- config$columns$description
ratio_col <- config$columns$gene_ratio
p_col <- config$columns$pvalue
n_col <- config$columns$count
group_col <- config$columns$group
has_group <- !is.null(group_col) && group_col %in% names(df)

df[[desc_col]] <- as.character(df[[desc_col]])
df[[ratio_col]] <- as.numeric(df[[ratio_col]])
df[[p_col]] <- as.numeric(df[[p_col]])
df[[n_col]] <- as.numeric(df[[n_col]])
if (any(!nzchar(df[[desc_col]]))) {
    stop("Description must be non-empty.", call. = FALSE)
}
ok <- is.finite(df[[ratio_col]]) & is.finite(df[[p_col]]) &
    df[[p_col]] > 0 & is.finite(df[[n_col]]) & df[[n_col]] > 0
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with non-finite ratio, p, or count.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

df$neg_log10 <- -log10(df[[p_col]])
if (identical(config$sort, "pvalue")) {
    df <- df[order(df$neg_log10, df[[desc_col]]), , drop = FALSE]
} else {
    df <- df[order(df[[ratio_col]], df[[desc_col]]), , drop = FALSE]
}
df[[desc_col]] <- factor(df[[desc_col]], levels = unique(df[[desc_col]]))
if (has_group) {
    df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
}

ink <- palette_colors("Qualitative.Safe")[[5]]
fill_cols <- palette_colors(config$palettes$fill)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[ratio_col]],
    y = .data[[desc_col]],
    fill = neg_log10,
    size = .data[[n_col]]
)) +
    geom_point(shape = 21, colour = ink, stroke = 0.35) +
    scale_fill_gradientn(
        name = config$labels$fill,
        colours = fill_cols
    ) +
    scale_size_continuous(
        name = config$labels$size,
        range = config$size_range
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    labs(x = config$labels$x, y = NULL) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        plot.margin = margin(4, 10, 4, 4, "mm"),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        axis.title.x = element_text(size = 11),
        axis.text = element_text(size = 9)
    )

if (isTRUE(config$facet) && has_group) {
    p <- p + facet_wrap(vars(.data[[group_col]]), scales = "free_y")
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p, io$output,
    width = config$size$width,
    height = config$size$height
)
