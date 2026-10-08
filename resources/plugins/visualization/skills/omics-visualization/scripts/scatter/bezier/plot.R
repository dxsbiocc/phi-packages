#!/usr/bin/env Rscript

# Template-ID: scatter-bezier
#
# Purpose:
#   Draw a paired scatter plot with Bezier curves joining observations
#   across two sample types, plus a paired Wilcoxon test.
#
# Inputs:
#   A table with one row per pair-by-sample-type. Default example:
#     - sample_type: sample label
#     - BRCA1: numeric expression
#     - patient: pairing identifier
#
# Output:
#   A PDF, PNG, or SVG Bezier paired scatter plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggforce, ggsignif, ggpubr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change sample order or how Bezier
#   control points are built from dodged point positions.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one paired observation at one sample type.
#   Bezier curves join the same pairing identifier across sample types.
#   Significance brackets use a paired Wilcoxon test at p < 0.05.

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
        x = "sample_type",
        y = "BRCA1",
        group = "patient"
    ),
    labels = list(
        title = "TCGA-ACC",
        x = "",
        y = "Expression of BRCA1"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggforce", "ggsignif", "ggpubr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group

df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

sig_comparisons <- list()
if (nlevels(df[[x_col]]) > 1) {
    stats_df <- compare_means(
        as.formula(paste(y_col, x_col, sep = " ~ ")),
        data = df,
        method = "wilcox.test",
        paired = TRUE
    )
    stats_df <- stats_df[stats_df$p < 0.05, , drop = FALSE]
    if (nrow(stats_df) > 0) {
        sig_comparisons <- Map(c, stats_df$group1, stats_df$group2)
    }
}

p_tmp <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]]
)) +
    geom_point(position = position_dodge2(width = 0.3))
built <- ggplot_build(p_tmp)$data[[1]]
n_first <- sum(built$group == 1)
built$index <- rep(seq_len(n_first), length(unique(built$group)))
bezier_data <- do.call(rbind, lapply(
    split(built, built$index),
    function(sub) {
        sub <- sub[order(sub$group), ]
        do.call(rbind, lapply(seq_len(nrow(sub) - 1), function(i) {
            data.frame(
                x = c(
                    sub$x[[i]], sub$x[[i]] + 0.3,
                    sub$x[[i + 1]] - 0.3, sub$x[[i + 1]]
                ),
                y = c(
                    sub$y[[i]], sub$y[[i]],
                    sub$y[[i + 1]], sub$y[[i + 1]]
                ),
                group = paste0(sub$index[[1]], i)
            )
        }))
    }
))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#2dde98", "#ff4f81")

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]]
)) +
    geom_blank() +
    geom_bezier(
        data = bezier_data,
        aes(x = x, y = y, group = group),
        colour = "#caccd1",
        alpha = 0.5,
        inherit.aes = FALSE
    ) +
    geom_point(
        position = position_dodge2(width = 0.3),
        size = 4,
        colour = "white"
    ) +
    geom_point(
        aes(fill = .data[[config$columns$x]]),
        position = position_dodge2(width = 0.3),
        size = 3,
        shape = 21,
        stroke = 0.5,
        alpha = 0.9,
        show.legend = FALSE
    )

if (length(sig_comparisons) > 0) {
    p <- p + geom_signif(
        comparisons = sig_comparisons,
        test = "wilcox.test",
        test.args = list(paired = TRUE),
        margin_top = 0.1,
        map_signif_level = FALSE
    )
}

p <- p +
    scale_fill_manual(values = fill_values) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
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
