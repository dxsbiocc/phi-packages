#!/usr/bin/env Rscript

# Template-ID: sankey-parallel-sets
#
# Purpose:
#   Draw a parallel-sets (alluvial) diagram: conserved counts flow
#   across several categorical axes.
#
# Inputs:
#   One row per combination of axis categories. Default example:
#     - tissue, cluster, response: axis categories
#     - n: supplied count
#
# Output:
#   A PDF, PNG, or SVG parallel-sets diagram.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggforce
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (axis columns, value, fill).
#   Edit DATA PREPARATION to change axis order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one already-counted combination of categories.
#   This script does not tabulate raw samples and does not test
#   association. Ribbon width is the supplied count.
#   Axis order is CONFIG, not a clustering of the categories.

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
        axes = c("tissue", "cluster", "response"),
        value = "n",
        fill = "cluster"
    ),
    axis_width = 0.12,
    alpha = 0.55,
    palette = "Brand.Algolia",
    labels = list(
        title = "",
        x = "",
        y = "",
        fill = NULL
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggforce"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
axis_cols <- config$columns$axes
val_col <- config$columns$value
fill_col <- config$columns$fill
require_columns(df, c(axis_cols, val_col, fill_col))

if (length(axis_cols) < 2L) {
    stop("Need at least two axis columns.", call. = FALSE)
}

df[[val_col]] <- as.numeric(df[[val_col]])
ok <- is.finite(df[[val_col]]) & df[[val_col]] > 0
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with non-positive or non-finite n.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

for (col in axis_cols) {
    df[[col]] <- as.character(df[[col]])
}
df[[fill_col]] <- factor(df[[fill_col]], levels = unique(df[[fill_col]]))

keep <- unique(c(axis_cols, val_col, fill_col))
df <- df[keep]
gathered <- gather_set_data(df, seq_along(axis_cols))
gathered$axis <- factor(
    gathered$x,
    levels = seq_along(axis_cols),
    labels = axis_cols
)

fill_n <- nlevels(gathered[[fill_col]])
fill_cols <- palette_colors(config$palette, n = fill_n)
# Light axis bars (RdBu midpoint) so Bold ribbons stay the focus.
axis_fill <- palette_colors("Diverging.RdBu")[[5]]
label_col <- palette_colors("Qualitative.Bold")[[11]]
fill_name <- config$labels$fill
if (is.null(fill_name)) {
    fill_name <- fill_col
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(gathered, aes(
    x = axis,
    id = id,
    split = y,
    value = .data[[val_col]]
)) +
    geom_parallel_sets(
        aes(fill = .data[[fill_col]]),
        alpha = config$alpha,
        axis.width = config$axis_width,
        sep = 0.05,
        strength = 0.5
    ) +
    geom_parallel_sets_axes(
        axis.width = config$axis_width,
        fill = axis_fill,
        colour = label_col,
        linewidth = 0.35
    ) +
    geom_parallel_sets_labels(
        colour = label_col,
        angle = 0,
        size = 3.2
    ) +
    scale_fill_manual(values = fill_cols, name = fill_name) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 9, height = 5.5)
