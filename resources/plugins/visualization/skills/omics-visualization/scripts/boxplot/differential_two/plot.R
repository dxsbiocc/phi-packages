#!/usr/bin/env Rscript

# Template-ID: boxplot-differential-two
#
# Purpose:
#   Draw grouped boxplots of gene expression for two groups, with
#   significance labels on each gene.
#
# Inputs:
#   A table with one row per observation. Default example:
#     - gene: category on the x-axis
#     - exprs: numeric expression
#     - group: fill group (two levels)
#
# Output:
#   A PDF, PNG, or SVG two-group differential boxplot.
#
# Dependencies:
#   ggplot2, readr, ggpubr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change gene or group order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one unpaired observation of a gene in a group.
#   Significance labels compare the two groups within each gene.
#   Non-significant comparisons are omitted.

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
        x = "gene",
        y = "exprs",
        fill = "group"
    ),
    labels = list(
        title = "",
        x = "",
        y = "Gene expression level"
    )
)

load_packages(c("ggplot2", "readr", "ggpubr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
fill_col <- config$columns$fill
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[fill_col]] <- factor(df[[fill_col]], levels = unique(df[[fill_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c(
    "#d7689d", "#f3d337", "#a0c443", "#66b5ae", "#5da4dc", "#40769e"
)

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[fill_col]]
)) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(
        aes(group = .data[[fill_col]]),
        label = "p.signif",
        hide.ns = TRUE
    ) +
    scale_fill_manual(values = fill_values) +
    guides(x = guide_prism_bracket()) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
