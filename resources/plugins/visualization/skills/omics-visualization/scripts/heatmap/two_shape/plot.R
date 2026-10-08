#!/usr/bin/env Rscript

# Template-ID: heatmap-two-shape
#
# Purpose:
#   Draw a split-shape correlation heatmap: rings for one group's
#   correlations on the lower triangle and stars for the other group
#   on the upper triangle.
#
# Inputs:
#   A table with numeric variable columns plus a two-level group
#   column. Default example:
#     - group: condition (two levels)
#     - remaining columns: variables correlated within each group
#
# Output:
#   A PDF, PNG, or SVG two-group shape correlation heatmap.
#
# Dependencies:
#   ggplot2, readr, ggcor
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change which group is lower vs upper
#   or the correlation method.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each group is correlated separately (Pearson). Lower-triangle
#   rings are the first group in sorted unique order; upper-triangle
#   stars are the second group. This is a display split, not a test
#   of group difference.

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
        group = "group"
    ),
    labels = list(
        title = "",
        x = "",
        y = "",
        fill = "Correlation"
    )
)

load_packages(c("ggplot2", "readr", "ggcor"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

group_col <- config$columns$group
feature_cols <- setdiff(names(df), group_col)
if (length(feature_cols) < 2L) {
    stop("Need at least two numeric columns besides the group column.",
        call. = FALSE
    )
}

groups <- sort(unique(as.character(df[[group_col]])))
if (length(groups) != 2L) {
    stop("group column must have exactly two levels.", call. = FALSE)
}

for (col in feature_cols) {
    df[[col]] <- as.numeric(df[[col]])
}

g1 <- as.data.frame(df[df[[group_col]] == groups[[1]], feature_cols])
g2 <- as.data.frame(df[df[[group_col]] == groups[[2]], feature_cols])

cor_tbl <- as_cor_tbl(
    correlate(g1, method = "pearson", cor.test = TRUE),
    extra.mat = list(r2 = stats::cor(g2, use = "pairwise.complete.obs"))
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_colours <- grDevices::colorRampPalette(c(
    "#ed6ca4", "#ffffffff", "#7bc4e2"
))(50)

p <- quickcor(cor_tbl) +
    geom_ring(
        mapping = aes(fill = r),
        data = get_data(type = "lower")
    ) +
    geom_star(
        mapping = aes(fill = r2),
        data = get_data(type = "upper")
    ) +
    geom_diag_label(size = 2.5) +
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
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        panel.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
