#!/usr/bin/env Rscript

# Template-ID: heatmap-shape
#
# Purpose:
#   Draw a lower-triangle correlation heatmap using star-shaped cells
#   instead of tiles.
#
# Inputs:
#   A numeric table (tsv/csv), one row per sample and one column per
#   variable. Default example uses ferroptosis-related genes.
#
# Output:
#   A PDF, PNG, or SVG shape-overlay correlation heatmap.
#
# Dependencies:
#   ggplot2, readr, linkET
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change the correlation method or which
#   columns are included.
#   Edit PLOT only when the marker shape or styling must change.
#
# Scientific assumptions:
#   Pairwise Pearson correlations are computed from paired rows.
#   Cell shape encodes the same correlation as fill; no p-values are
#   shown unless PLOT is changed.

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
    columns = list(),
    labels = list(
        title = "",
        x = "",
        y = "",
        fill = "Correlation"
    )
)

load_packages(c("ggplot2", "readr", "linkET"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

feature_cols <- names(df)
if (length(feature_cols) < 2L) {
    stop("Need at least two numeric columns to compute correlations.",
        call. = FALSE
    )
}
for (col in feature_cols) {
    df[[col]] <- as.numeric(df[[col]])
}

cor_obj <- correlate(
    as.data.frame(df[feature_cols]),
    method = "pearson"
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_colours <- grDevices::colorRampPalette(c(
    "#acd372", "#ffffffff", "#ed6ca4"
))(50)

p <- qcorrplot(cor_obj, type = "lower") +
    geom_shaping(marker = marker("star")) +
    scale_fill_gradientn(
        name = config$labels$fill,
        colours = fill_colours
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme(
        plot.background = element_blank(),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_blank(),
        legend.position = "right",
        panel.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
