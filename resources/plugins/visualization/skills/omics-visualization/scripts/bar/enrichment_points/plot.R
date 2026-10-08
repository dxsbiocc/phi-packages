#!/usr/bin/env Rscript

# Template-ID: bar-enrichment-points
#
# Purpose:
#   Draw a faceted enrichment bar chart with gene-ratio points overlaid
#   on a shared numeric axis.
#
# Inputs:
#   A table with one row per term. Default example:
#     - Description: term name
#     - ONTOLOGY: term category
#     - GeneRatio: gene ratio
#     - qvalue: adjusted q-value
#
# Output:
#   A PDF, PNG, or SVG enrichment bar-and-point chart.
#
# Dependencies:
#   ggplot2, readr, ggprism, gground
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change sort order or the -log10 / percent
#   transforms.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one enrichment term with a q-value and gene ratio.
#   -log10(qvalue) and GeneRatio*100 are display transforms.
#   Terms are sorted within each category by q-value.

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
        group = "ONTOLOGY",
        gene_ratio = "GeneRatio",
        qvalue = "qvalue"
    ),
    labels = list(
        title = "",
        x = "Gene Ratio(%) \t -log10(p.adjust)",
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
group_col <- config$columns$group
ratio_col <- config$columns$gene_ratio
q_col <- config$columns$qvalue

df[[q_col]] <- as.numeric(df[[q_col]])
df[[ratio_col]] <- as.numeric(df[[ratio_col]])
df[[desc_col]] <- as.character(df[[desc_col]])
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
min_pos <- min(df[[q_col]][df[[q_col]] > 0], na.rm = TRUE)
df$neg_log10 <- -log10(pmax(df[[q_col]], min_pos / 10))
df$gene_pct <- -df[[ratio_col]] * 100

df <- df[order(df[[group_col]], df[[q_col]]), , drop = FALSE]
df[[desc_col]] <- factor(df[[desc_col]], levels = unique(df[[desc_col]]))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#3BE8B0", "#1AAFD0", "#6A67CE", "#FFB900", "#FC636B")
colour_values <- c("#3be8b0", "#1aafd0", "#6a67ce", "#ffb900", "#fc636b")

p <- ggplot(df, aes(
    x = neg_log10,
    y = .data[[desc_col]],
    fill = .data[[group_col]]
)) +
    geom_round_col(alpha = 0.6, show.legend = FALSE) +
    geom_line(
        aes(x = gene_pct, group = 1),
        orientation = "y"
    ) +
    geom_point(
        aes(x = gene_pct),
        colour = "white",
        size = 4,
        show.legend = FALSE
    ) +
    geom_point(
        aes(x = gene_pct, colour = .data[[group_col]]),
        size = 3,
        show.legend = FALSE
    ) +
    geom_text(
        aes(x = 0, label = .data[[desc_col]]),
        hjust = 0,
        size = 3,
        lineheight = 0.7
    ) +
    scale_fill_manual(name = "Category", values = fill_values) +
    scale_colour_manual(values = colour_values) +
    scale_x_continuous(limits = c(-9, 35), expand = c(0, 0)) +
    facet_wrap(vars(.data[[group_col]]), scales = "free") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        legend.title = element_blank(),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
