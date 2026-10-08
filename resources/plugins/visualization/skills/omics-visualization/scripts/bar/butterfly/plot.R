#!/usr/bin/env Rscript

# Template-ID: bar-butterfly
#
# Purpose:
#   Draw a butterfly (back-to-back) bar chart comparing a numeric value
#   across two groups on a shared categorical axis.
#
# Inputs:
#   A table with one row per category-by-group. Default example:
#     - Age: category label
#     - Sex: group (two levels)
#     - Count: numeric value
#
# Output:
#   A PDF, PNG, or SVG butterfly bar chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to choose which group goes left or to change
#   category order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category-by-group pair with a comparable numeric value.
#   Values are treated as already summarized, not as raw replicates.
#   Negating one group's values is a display transform, not a data transform.

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
        x = "Count",
        y = "Age",
        group = "Sex"
    ),
    labels = list(
        title = "",
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
group_col <- config$columns$group

df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df[[x_col]] <- as.numeric(df[[x_col]])

# First group in file order is drawn to the left.
left_group <- levels(df[[group_col]])[[1]]
df[[x_col]] <- ifelse(
    df[[group_col]] == left_group,
    -df[[x_col]],
    df[[x_col]]
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c(
    "#d7689d", "#f3d337", "#a0c443", "#66b5ae", "#5da4dc", "#40769e"
)

x_max <- max(abs(df[[x_col]]), na.rm = TRUE)
x_breaks <- pretty(c(-x_max, x_max), n = 5)

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[group_col]]
)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(
        vars(.data[[group_col]]),
        scales = "free_x",
        nrow = 1,
        axes = "all_y"
    ) +
    scale_fill_manual(values = fill_values) +
    scale_x_continuous(
        breaks = x_breaks,
        labels = abs(x_breaks)
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12, hjust = 0.5),
        panel.spacing = unit(0, "lines"),
        plot.background = element_blank()
    )

gt <- ggplotGrob(p)
idx <- which(gt$layout$name == "axis-l-1-1")
if (length(idx) == 1L) {
    left <- gt$layout$l[[idx]]
    gt$grobs[[idx]] <- grid::nullGrob()
    gt$widths[left] <- unit(0, "cm")
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(gt, io$output)
