#!/usr/bin/env Rscript

# Template-ID: pie-nest
#
# Purpose:
#   Draw a nested (two-ring) pie chart with inner and outer series.
#
# Inputs:
#   A table with one row per slice. Default example:
#     - name: slice label
#     - value: numeric size
#     - category: ring (inner / outer)
#
# Output:
#   A PDF, PNG, or SVG nested pie chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or reshape the table.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Rows tagged inner form the inner ring; outer form the outer ring.
#   The two rings are separate series, not a parent-child rollup computed here.
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
        x = "name",
        y = "value",
        ring = "category"
    ),
    labels = list(
        title = "Nested Pie Chart",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
ring_col <- config$columns$ring
df[[y_col]] <- as.numeric(df[[y_col]])
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df$ring_x <- ifelse(df[[ring_col]] == unique(df[[ring_col]])[[1]], 1, 2)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = ring_x,
    y = .data[[y_col]],
    fill = .data[[x_col]]
)) +
    geom_col(width = 1, colour = "white") +
    scale_fill_manual(
        values = expand_palette(
            c("#f79a3e",
            "#00363d",
            "#eb4962",
            "#37b8af",
            "#78a300",
            "#f0ca28",
            "#30aabc",
            "#eb6651",
            "#4c5058",
            "#6e7079"),
            nlevels(df[[x_col]])
        )
    ) +
    coord_polar(theta = "y") +
    xlim(0, 2.6) +
    labs(title = config$labels$title, x = NULL, y = NULL, fill = "") +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 8)
