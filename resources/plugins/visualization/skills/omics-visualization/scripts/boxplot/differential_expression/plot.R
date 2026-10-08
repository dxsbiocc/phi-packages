#!/usr/bin/env Rscript

# Template-ID: boxplot-differential-expression
#
# Purpose:
#   Draw a boxplot of a numeric value across categories, with jittered
#   points and Wilcoxon significance brackets for pairs at p < 0.05.
#
# Inputs:
#   A table with one row per sample. Default example:
#     - stage: category label
#     - TP53: numeric value
#
# Output:
#   A PDF, PNG, or SVG differential boxplot.
#
# Dependencies:
#   ggplot2, readr, ggsignif, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or which pairs
#   are tested.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one unpaired observation in a category.
#   Brackets are Wilcoxon rank-sum tests at p < 0.05, without
#   multiple-testing correction.
#   Non-significant pairs are omitted.

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
        x = "stage",
        y = "TP53"
    ),
    labels = list(
        title = "TCGA-ACC",
        x = "",
        y = "Expression of TP53"
    )
)

load_packages(c("ggplot2", "readr", "ggsignif", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

x_lvls <- levels(df[[x_col]])
comparisons <- list()
for (i in seq_len(length(x_lvls) - 1L)) {
    for (j in seq(i + 1L, length(x_lvls))) {
        a <- df[[y_col]][df[[x_col]] == x_lvls[[i]]]
        b <- df[[y_col]][df[[x_col]] == x_lvls[[j]]]
        pv <- tryCatch(
            wilcox.test(a, b, paired = FALSE)$p.value,
            error = function(e) 1
        )
        if (is.finite(pv) && pv < 0.05) {
            comparisons[[length(comparisons) + 1L]] <- c(
                x_lvls[[i]], x_lvls[[j]]
            )
        }
    }
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#337ab7", "#5bc0de", "#5cb85c", "#f0ad4e", "#d9534f")

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[x_col]]
)) +
    geom_boxplot(
        alpha = 0.8,
        outlier.shape = NA,
        show.legend = FALSE
    ) +
    geom_point(
        position = position_jitter(width = 0.2),
        size = 3,
        alpha = 0.6,
        stroke = 0,
        show.legend = FALSE
    ) +
    scale_fill_manual(values = fill_values) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
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

if (length(comparisons) > 0) {
    p <- p + geom_signif(
        comparisons = comparisons,
        test = "wilcox.test",
        test.args = list(paired = FALSE),
        margin_top = 0.1,
        step_increase = 0.08,
        map_signif_level = FALSE
    )
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
