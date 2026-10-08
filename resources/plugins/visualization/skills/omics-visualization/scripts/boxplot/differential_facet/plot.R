#!/usr/bin/env Rscript

# Template-ID: boxplot-differential-facet
#
# Purpose:
#   Draw faceted boxplots of gene expression across two groups, with
#   jittered points and significance labels in each panel.
#
# Inputs:
#   A table with one row per observation. Default example:
#     - group: comparison group on the x-axis
#     - exprs: numeric expression
#     - gene: facet panel
#
# Output:
#   A PDF, PNG, or SVG faceted differential boxplot.
#
# Dependencies:
#   ggplot2, readr, ggpubr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change group or facet order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one unpaired observation of a gene in a group.
#   Facets are independent panels; tests compare groups within a gene.
#   Significance labels are pairwise tests with non-significant
#   pairs omitted.

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
        x = "group",
        y = "exprs",
        colour = "group",
        facet = "gene"
    ),
    labels = list(
        title = "",
        x = "",
        y = "Gene expression level"
    )
)

load_packages(c("ggplot2", "readr", "ggpubr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
colour_col <- config$columns$colour
facet_col <- config$columns$facet
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[colour_col]] <- factor(df[[colour_col]], levels = unique(df[[colour_col]]))
df[[facet_col]] <- factor(df[[facet_col]], levels = unique(df[[facet_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c(
    "#8F499C", "#4185BE", "#6DC067", "#F6DB35", "#F78822"
)

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    colour = .data[[colour_col]]
)) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    geom_point(
        position = position_jitterdodge(
            jitter.width = 0.15,
            dodge.width = 0.7
        ),
        size = 2,
        alpha = 0.6,
        show.legend = FALSE
    ) +
    stat_compare_means(
        aes(group = .data[[x_col]]),
        label = "p.signif",
        hide.ns = TRUE,
        show.legend = FALSE
    ) +
    scale_colour_manual(values = colour_values) +
    facet_wrap(
        vars(.data[[facet_col]]),
        nrow = 1,
        scales = "free_x"
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_bw() +
    theme(
        plot.background = element_blank(),
        axis.text = element_text(face = "bold", size = 8),
        axis.title.y = element_text(size = 10, face = "bold"),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "#6dc067"),
        panel.grid = element_blank(),
        panel.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
