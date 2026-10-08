#!/usr/bin/env Rscript

# Template-ID: bar-polar-horizontal
#
# Purpose:
#   Draw a polar bar chart with values mapped to angle (a circular stacked ring).
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - day: category label
#     - direct: numeric value
#
# Output:
#   A PDF, PNG, or SVG polar (angular) bar chart.
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
#   Each row is one category with a comparable numeric value.
#   Mapping value to angle is a display choice, not a statistical transform.
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
        x = "day",
        y = "direct"
    ),
    labels = list(
        title = "Polar Horizontal Bar Chart",
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
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])
df <- df[order(df[[x_col]]), , drop = FALSE]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = 1,
    y = .data[[y_col]],
    fill = .data[[x_col]]
)) +
    geom_col(width = 0.5, colour = "white") +
    scale_fill_manual(
        values = expand_palette(
            c("#5070dd",
            "#b6d634",
            "#505372",
            "#ff994d",
            "#0ca8df",
            "#ffd10a",
            "#fb628b"),
            nlevels(df[[x_col]])
        )
    ) +
    coord_polar(theta = "y") +
    xlim(0.4, 1.5) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        fill = ""
    ) +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 8)
