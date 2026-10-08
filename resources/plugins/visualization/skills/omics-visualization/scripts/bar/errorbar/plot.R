#!/usr/bin/env Rscript

# Template-ID: bar-errorbar
#
# Purpose:
#   Draw grouped bars of mean values with standard-error whiskers and
#   jittered raw points.
#
# Inputs:
#   A table with one row per observation. Default example:
#     - condition: x-axis group
#     - value: numeric measurement
#     - specie: fill group
#
# Output:
#   A PDF, PNG, or SVG bar chart with error bars.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change group order.
#   Edit PLOT to change the summary function or error-bar statistic.
#
# Scientific assumptions:
#   Each row is one replicate observation, not a pre-aggregated mean.
#   Bars are the mean; whiskers are mean ± standard error.

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
        x = "condition",
        y = "value",
        fill = "specie"
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
fill_col <- config$columns$fill
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[fill_col]] <- factor(df[[fill_col]], levels = unique(df[[fill_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

mean_se <- function(x) {
    x <- x[is.finite(x)]
    m <- mean(x)
    se <- stats::sd(x) / sqrt(length(x))
    data.frame(y = m, ymin = m - se, ymax = m + se)
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c(
    "#009392", "#39b185", "#9ccb86", "#e9e29c",
    "#eeb479", "#e88471", "#cf597e"
)

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[fill_col]]
)) +
    stat_summary(
        fun = mean,
        geom = "col",
        position = position_dodge(width = 0.9),
        width = 0.7,
        colour = "grey30"
    ) +
    stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        position = position_dodge(width = 0.9),
        width = 0.2
    ) +
    geom_jitter(
        position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
        alpha = 0.7,
        size = 3
    ) +
    scale_fill_manual(values = fill_values) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(plot.background = element_blank())

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
