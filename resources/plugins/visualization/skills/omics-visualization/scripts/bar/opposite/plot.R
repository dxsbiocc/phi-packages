#!/usr/bin/env Rscript

# Template-ID: bar-opposite
#
# Purpose:
#   Draw a two-group reverse bar chart of enrichment terms, with one
#   group to the right and the other to the left of a shared axis.
#
# Inputs:
#   A table with one row per term. Default example:
#     - Description: term name
#     - p.adjust: adjusted p-value
#     - Change: group (exactly two levels)
#
# Output:
#   A PDF, PNG, or SVG opposite bar chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to choose which group goes right or to change
#   sort order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one enrichment term with an adjusted p-value.
#   Exactly two groups are required. The first group in file order is
#   drawn to the right.
#   -log10(p.adjust) is a display transform; the sign flip is visual.

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
        description = "Description",
        pvalue = "p.adjust",
        group = "Change"
    ),
    labels = list(
        title = "",
        x = "-log10(p.adjust)",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

desc_col <- config$columns$description
p_col <- config$columns$pvalue
group_col <- config$columns$group

df[[p_col]] <- as.numeric(df[[p_col]])
df[[desc_col]] <- as.character(df[[desc_col]])
df[[group_col]] <- as.character(df[[group_col]])
groups <- unique(df[[group_col]])
if (length(groups) != 2) {
    stop("Only support two groups!", call. = FALSE)
}

df[[group_col]] <- factor(df[[group_col]], levels = groups)
min_pos <- min(df[[p_col]][df[[p_col]] > 0], na.rm = TRUE)
df$neg_log10 <- -log10(pmax(df[[p_col]], min_pos / 10))
df$sign <- ifelse(df[[group_col]] == groups[[1]], 1, -1)
df$margin <- ifelse(df$sign == 1, -0.2, 0.2)
df$hjust <- ifelse(df$sign == 1, 1, 0)
df$signed_x <- df$neg_log10 * df$sign

df <- df[order(df[[group_col]], df$neg_log10), , drop = FALSE]
df[[desc_col]] <- factor(df[[desc_col]], levels = unique(df[[desc_col]]))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#acd372", "#ed6ca4")

p <- ggplot(df, aes(
    x = signed_x,
    y = .data[[desc_col]],
    fill = .data[[group_col]]
)) +
    geom_col(alpha = 0.7) +
    geom_text(
        aes(x = margin, label = .data[[desc_col]], hjust = hjust),
        colour = "black",
        show.legend = FALSE
    ) +
    scale_fill_manual(name = "Correlation", values = fill_values) +
    scale_x_continuous(
        limits = c(-19, 25),
        breaks = c(-20, -10, 0, 10, 20),
        labels = c(20, 10, 0, 10, 20)
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top",
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
