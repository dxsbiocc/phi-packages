#!/usr/bin/env Rscript

# Template-ID: bar-percent
#
# Purpose:
#   Draw faceted horizontal percent bars comparing a numeric frequency
#   across categories, one panel per group.
#
# Inputs:
#   A table with one row per category-by-group. Default example:
#     - type: category label
#     - symbol: group (e.g. gene)
#     - freq: numeric frequency
#
# Output:
#   A PDF, PNG, or SVG grouped percent bar chart.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggh4x
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change the percent transform or category
#   order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category-by-group pair with a comparable frequency.
#   Negating freq*100 is a display transform so bars grow leftward.
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
        y = "type",
        x = "freq",
        group = "symbol"
    ),
    labels = list(
        title = "",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggh4x"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group

df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df[[x_col]] <- as.numeric(df[[x_col]])
df$pct <- -df[[x_col]] * 100

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c(
    "#009392", "#39b185", "#9ccb86", "#e9e29c", "#eeb479", "#e88471", "#cf597e"
)
n_group <- nlevels(df[[group_col]])
fill_values <- rep_len(fill_values, n_group)

p <- ggplot(df, aes(
    x = pct,
    y = .data[[y_col]],
    fill = .data[[group_col]]
)) +
    geom_col(aes(x = -0.5), alpha = 0.2) +
    geom_col() +
    facet_wrap2(
        vars(.data[[group_col]]),
        nrow = 1,
        strip = strip_themed(
            text_x = elem_list_text(size = 15, vjust = 0),
            background_x = elem_list_rect(fill = fill_values)
        )
    ) +
    scale_fill_manual(name = "Genes", values = fill_values) +
    scale_x_continuous(
        expand = c(0, 0),
        limits = c(-0.5, 0),
        breaks = c(-0.5, 0),
        labels = c("50%", "0%")
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    guides(fill = "none") +
    theme_prism() +
    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.background = element_blank()
    )

gt <- ggplotGrob(p)
axis_b_indices <- grep("axis-b", gt$layout$name)
if (length(axis_b_indices) > 1) {
    for (i in axis_b_indices[-1]) {
        gt$grobs[[i]] <- grid::nullGrob()
    }
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(gt, io$output)
