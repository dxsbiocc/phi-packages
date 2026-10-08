#!/usr/bin/env Rscript

# Template-ID: scatter-marginal
#
# Purpose:
#   Draw a grouped scatter plot with a linear fit, Pearson correlation,
#   and boxplot margins on both axes.
#
# Inputs:
#   A table with one row per observation. Default example:
#     - x: numeric x value
#     - y: numeric y value
#     - group: grouping label
#
# Output:
#   A PDF, PNG, or SVG scatter plot with marginal boxplots.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggExtra, ggpubr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change group order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one observation with two numeric variables.
#   The linear fit and Pearson correlation are descriptive overlays.
#   Marginal boxplots summarize the same points, not a separate sample.

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
        x = "x",
        y = "y",
        group = "group"
    ),
    labels = list(
        title = "",
        x = "X value",
        y = "Y value"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggExtra", "ggpubr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group

df <- df[stats::complete.cases(df[, c(x_col, y_col, group_col)]), ]
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#333463", "#079C9C")

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]],
    colour = .data[[config$columns$group]],
    fill = .data[[config$columns$group]]
)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.8) +
    stat_cor(method = "pearson") +
    scale_fill_manual(values = fill_values) +
    scale_colour_manual(values = fill_values) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none"
    )

p <- suppressWarnings(ggMarginal(
    p,
    type = "boxplot",
    groupFill = TRUE,
    margins = "both",
    groupColour = TRUE
))

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
