#!/usr/bin/env Rscript

# Template-ID: scatter-correlation
#
# Purpose:
#   Draw a scatter plot of two numeric variables with a linear fit
#   and a Pearson correlation label.
#
# Inputs:
#   A table with one row per observation. Default example:
#     - BRCA1: first numeric axis
#     - BRCA2: second numeric axis
#
# Output:
#   A PDF, PNG, or SVG correlation scatter plot.
#
# Dependencies:
#   ggplot2, readr, ggpubr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change how the colour scale is derived.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one paired observation of two numeric variables.
#   The trend line is an ordinary least-squares fit.
#   The labelled statistic is a Pearson correlation, not a causal claim.

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
        x = "BRCA1",
        y = "BRCA2"
    ),
    labels = list(
        title = "",
        x = "BRCA1",
        y = "BRCA2"
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
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df$colour_sum <- df[[x_col]] + df[[y_col]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]]
)) +
    geom_point(size = 4.5, colour = "black") +
    geom_point(
        aes(colour = colour_sum),
        size = 4,
        show.legend = FALSE
    ) +
    geom_smooth(
        method = "lm",
        formula = y ~ x,
        colour = "#009392",
        fill = "#B1C7B3"
    ) +
    stat_cor(method = "pearson", size = 6) +
    scale_colour_gradient2(
        low = "#F1EAC8",
        mid = "#E5B9AD",
        high = "#D0587E"
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
