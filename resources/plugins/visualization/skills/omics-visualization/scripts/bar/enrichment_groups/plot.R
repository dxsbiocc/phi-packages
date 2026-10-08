#!/usr/bin/env Rscript

# Template-ID: bar-enrichment-groups
#
# Purpose:
#   Draw a multi-group enrichment bar chart with gene lists, count
#   bubbles, and category labels along the left side.
#
# Inputs:
#   A table with one row per term. Default example:
#     - index: row order
#     - ONTOLOGY: term category
#     - Description: term name
#     - p.adjust: adjusted p-value
#     - geneID: slash-separated gene symbols
#     - Count: gene count
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
#   Edit DATA PREPARATION to change category blocks or the -log10
#   transform.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one enrichment term with an adjusted p-value and count.
#   -log10(p.adjust) is a display transform, not an upstream test.
#   Left-side category rectangles are visual labels, not additional terms.

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
        y = "index",
        group = "ONTOLOGY",
        description = "Description",
        pvalue = "p.adjust",
        gene_id = "geneID",
        count = "Count"
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

y_col <- config$columns$y
group_col <- config$columns$group
desc_col <- config$columns$description
p_col <- config$columns$pvalue
gene_col <- config$columns$gene_id
count_col <- config$columns$count

df[[p_col]] <- as.numeric(df[[p_col]])
df[[count_col]] <- as.numeric(df[[count_col]])
df[[desc_col]] <- as.character(df[[desc_col]])
df[[gene_col]] <- as.character(df[[gene_col]])
df[[group_col]] <- as.character(df[[group_col]])
df <- df[order(df[[y_col]]), , drop = FALSE]
df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
min_pos <- min(df[[p_col]][df[[p_col]] > 0], na.rm = TRUE)
df$neg_log10 <- -log10(pmax(df[[p_col]], min_pos / 10))

point_x <- 0.5
width <- point_x
group_n <- as.integer(table(df[[group_col]])[levels(df[[group_col]])])
ymax <- cumsum(group_n)
ymin <- c(0, ymax[-length(ymax)]) + 0.6
ymax <- ymax + 0.4
rect_data <- data.frame(
    xmin = -3 * width,
    xmax = -2 * width,
    ymin = ymin,
    ymax = ymax,
    stringsAsFactors = FALSE
)
rect_data[[group_col]] <- levels(df[[group_col]])
rect_data$label_x <- (rect_data$xmin + rect_data$xmax) / 2
rect_data$label_y <- (rect_data$ymin + rect_data$ymax) / 2

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#7bc4e2", "#acd372", "#fbb05b", "#ed6ca4")

p <- ggplot(df, aes(
    x = neg_log10,
    y = .data[[y_col]],
    fill = .data[[group_col]]
)) +
    geom_round_col(width = 0.6, alpha = 0.8) +
    geom_text(
        aes(x = 0.05, label = .data[[desc_col]]),
        hjust = 0,
        size = 5
    ) +
    geom_text(
        aes(
            x = 0.1,
            label = .data[[gene_col]],
            colour = .data[[group_col]]
        ),
        hjust = 0,
        vjust = 2.6,
        size = 3.5,
        fontface = "italic",
        show.legend = FALSE
    ) +
    geom_point(
        aes(x = -point_x, size = .data[[count_col]]),
        shape = 21
    ) +
    geom_text(
        aes(x = -point_x, label = .data[[count_col]])
    ) +
    geom_round_rect(
        data = rect_data,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            fill = .data[[group_col]]
        ),
        inherit.aes = FALSE
    ) +
    geom_text(
        data = rect_data,
        aes(
            x = label_x,
            y = label_y,
            label = .data[[group_col]]
        ),
        inherit.aes = FALSE
    ) +
    annotate(
        "segment",
        x = 0, y = 0, xend = 8.6, yend = 0,
        linewidth = 1.5
    ) +
    scale_size_continuous(name = "Count", range = c(5, 16)) +
    scale_fill_manual(name = "Category", values = fill_values) +
    scale_colour_manual(values = fill_values) +
    scale_x_continuous(
        breaks = c(0, 2, 4, 6, 8),
        expand = expansion(mult = c(0, 0))
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 16, height = 10)
