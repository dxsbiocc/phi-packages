#!/usr/bin/env Rscript

# Template-ID: bar-radial
#
# Purpose:
#   Draw a radial (polar) bar chart of numeric values around a circle.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - type: category label
#     - n: numeric value
#
# Output:
#   A PDF, PNG, or SVG radial bar chart.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category with a comparable numeric value.
#   Polar coordinates are a display choice, not a statistical transform.
#   Values are treated as already summarized, not as raw replicates.

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
        x = "type",
        y = "n"
    ),
    labels = list(
        title = "",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c(
    "#d7689d", "#f3d337", "#a0c443", "#66b5ae", "#5da4dc", "#40769e"
)
n_fill <- nlevels(df[[x_col]])
if (length(fill_values) < n_fill) {
    fill_values <- grDevices::colorRampPalette(fill_values)(n_fill)
}

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[x_col]]
)) +
    geom_hline(yintercept = 0, color = "grey", linetype = "solid", linewidth = 1) +
    geom_hline(yintercept = 500, color = "grey", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = 1000, color = "grey", linetype = "dashed", linewidth = 1) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = fill_values) +
    scale_y_continuous(
        limits = c(-500, 1100),
        breaks = c(-500, 0, 500, 1000)
    ) +
    coord_polar() +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_minimal() +
    theme(
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
