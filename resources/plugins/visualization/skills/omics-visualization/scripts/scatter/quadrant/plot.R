#!/usr/bin/env Rscript

# Template-ID: scatter-quadrant
#
# Purpose:
#   Draw a quadrant scatter plot of two log2 fold-change axes, with
#   threshold bands and labels for the strongest points in each quadrant.
#
# Inputs:
#   A table with one row per gene. Default example:
#     - OE: overexpression log2 fold-change
#     - KO: knockout log2 fold-change
#     - group: quadrant label (None / L-H / H-L / H-H / L-L)
#     - symbol: gene label
#
# Output:
#   A PDF, PNG, or SVG quadrant scatter plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggrepel
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change the fold-change threshold or how
#   many labels are drawn per quadrant.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one gene with two comparable log2 fold-change values.
#   |log2FC| > 1 on each axis is a display threshold, not a statistical test.
#   Points labelled "None" are background and are not highlighted.

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
        x = "OE",
        y = "KO",
        group = "group",
        label = "symbol"
    ),
    labels = list(
        title = "Yap1 Overexpression and Knockout",
        x = "Overexpression log2FC",
        y = "Knockout log2FC"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggrepel"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group
label_col <- config$columns$label

df <- df[stats::complete.cases(df[, c(x_col, y_col, group_col, label_col)]), ]
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))

fc_threshold <- 1
rect_data <- data.frame(
    xmin = c(-Inf, fc_threshold, -Inf, fc_threshold),
    xmax = c(-fc_threshold, Inf, -fc_threshold, Inf),
    ymin = c(fc_threshold, fc_threshold, -Inf, -Inf),
    ymax = c(Inf, Inf, -fc_threshold, -fc_threshold),
    group = factor(
        c("L-H", "H-H", "L-L", "H-L"),
        levels = levels(df[[group_col]])
    )
)

use_data <- df[tolower(as.character(df[[group_col]])) != "none", , drop = FALSE]
text_data <- do.call(rbind, lapply(
    split(use_data, use_data[[group_col]], drop = TRUE),
    function(sub) {
        score <- abs(sub[[x_col]]) + abs(sub[[y_col]])
        sub[head(order(score, decreasing = TRUE), 5), , drop = FALSE]
    }
))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c(
    "L-H" = "#3288bd",
    "H-H" = "#66c2a5",
    "L-L" = "#f46d43",
    "H-L" = "#d53e4f"
)

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]]
)) +
    geom_point(size = 1, alpha = 0.5, colour = "grey50") +
    geom_vline(
        xintercept = c(-fc_threshold, fc_threshold),
        colour = "grey50",
        linetype = 2,
        linewidth = 0.6
    ) +
    geom_hline(
        yintercept = c(-fc_threshold, fc_threshold),
        colour = "grey50",
        linetype = 2,
        linewidth = 0.6
    ) +
    geom_rect(
        data = rect_data,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            fill = group
        ),
        alpha = 0.4,
        show.legend = FALSE,
        inherit.aes = FALSE
    ) +
    geom_point(
        data = use_data,
        aes(
            x = .data[[config$columns$x]],
            y = .data[[config$columns$y]],
            fill = .data[[config$columns$group]],
            size = abs(.data[[config$columns$x]]) + abs(.data[[config$columns$y]])
        ),
        shape = 21,
        stroke = 0,
        show.legend = FALSE
    ) +
    geom_text_repel(
        data = text_data,
        aes(label = .data[[config$columns$label]]),
        size = 4,
        max.overlaps = 1000
    ) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_size_continuous(range = c(2, 6)) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        legend.title = element_text(size = 12)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
