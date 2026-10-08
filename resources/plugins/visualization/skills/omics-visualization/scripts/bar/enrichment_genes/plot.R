#!/usr/bin/env Rscript

# Template-ID: bar-enrichment-genes
#
# Purpose:
#   Draw a grouped horizontal bar chart of enrichment terms, with gene
#   lists shown under each term.
#
# Inputs:
#   A table with one row per term. Default example:
#     - Description: term name
#     - p.adjust: adjusted p-value
#     - geneID: slash-separated gene symbols
#     - Change: term group (e.g. Positive / Negative)
#
# Output:
#   A PDF, PNG, or SVG enrichment bar chart.
#
# Dependencies:
#   ggplot2, readr, ggprism, gground
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change sort order or the -log10 transform.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one enrichment term with an adjusted p-value and gene list.
#   -log10(p.adjust) is a display transform, not an upstream test.
#   Terms are sorted within each group by -log10(p.adjust).

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
        gene_id = "geneID",
        group = "Change"
    ),
    labels = list(
        title = "",
        x = "-log10(p.adjust)",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "gground"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

desc_col <- config$columns$description
p_col <- config$columns$pvalue
gene_col <- config$columns$gene_id
group_col <- config$columns$group

df[[p_col]] <- as.numeric(df[[p_col]])
df[[desc_col]] <- as.character(df[[desc_col]])
df[[gene_col]] <- as.character(df[[gene_col]])
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
min_pos <- min(df[[p_col]][df[[p_col]] > 0], na.rm = TRUE)
df$neg_log10 <- -log10(pmax(df[[p_col]], min_pos / 10))

df <- df[order(df[[group_col]], df$neg_log10), , drop = FALSE]
df[[desc_col]] <- factor(df[[desc_col]], levels = unique(df[[desc_col]]))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#79ceb8", "#e95f5c")

p <- ggplot(df, aes(
    x = neg_log10,
    y = .data[[desc_col]],
    fill = .data[[group_col]]
)) +
    geom_round_col(width = 0.6) +
    geom_text(
        aes(x = 0.1, label = .data[[desc_col]]),
        hjust = 0,
        size = 4
    ) +
    geom_text(
        aes(
            x = 0.1,
            label = .data[[gene_col]],
            colour = .data[[group_col]]
        ),
        hjust = 0,
        vjust = 2.8,
        size = 3.5,
        fontface = "italic",
        show.legend = FALSE
    ) +
    scale_fill_manual(name = "Correlation", values = fill_values) +
    scale_colour_manual(name = "Correlation", values = fill_values) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0.1, 0.1)) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.direction = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.1),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
