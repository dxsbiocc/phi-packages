#!/usr/bin/env Rscript

# Template-ID: bar-enrichment-expand
#
# Purpose:
#   Draw a grouped horizontal bar chart of enrichment terms, with category
#   headers expanded in the y-axis so pathway hierarchy is readable.
#
# Inputs:
#   A table with one row per term. Default example:
#     - category: term group
#     - Description: term name
#     - p.adjust: adjusted p-value
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
#   Edit DATA PREPARATION to change group spacing, sort order, or
#   the -log10 transform.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one enrichment term with an adjusted p-value.
#   -log10(p.adjust) is a display transform, not an upstream test.
#   Spacer rows between groups are visual, not additional terms.

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
        group = "category"
    ),
    labels = list(
        title = "KEGG enrichment",
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
group_col <- config$columns$group

df[[p_col]] <- as.numeric(df[[p_col]])
df[[group_col]] <- as.character(df[[group_col]])
df[[desc_col]] <- as.character(df[[desc_col]])
groups <- unique(df[[group_col]])

header <- df[rep(NA_integer_, length(groups)), , drop = FALSE]
header[[desc_col]] <- groups
header[[group_col]] <- groups
header[[p_col]] <- NA_real_

df <- rbind(header, df)
df[[group_col]] <- factor(df[[group_col]], levels = groups)
df <- df[order(
    df[[group_col]], !is.na(df[[p_col]]),
    df[[p_col]]
), , drop = FALSE]
df[[desc_col]] <- factor(df[[desc_col]], levels = unique(df[[desc_col]]))
df$neg_log10 <- -log10(pmax(df[[p_col]], 1e-300))
df$text_size <- ifelse(is.na(df[[p_col]]), 4, 3)
df$padding <- ifelse(is.na(df[[p_col]]), 0, 0.1)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#3BE8B0", "#1AAFD0", "#6A67CE", "#FFB900", "#FC636B")

p <- ggplot(df, aes(
    x = neg_log10,
    y = .data[[desc_col]],
    fill = .data[[group_col]]
)) +
    geom_round_col(orientation = "y", alpha = 0.6) +
    geom_text(
        aes(x = padding, label = .data[[desc_col]], size = text_size),
        hjust = 0,
        show.legend = FALSE
    ) +
    scale_size_continuous(range = c(4, 6)) +
    scale_fill_manual(name = "Category", values = fill_values) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    guides(size = "none", fill = "none") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(size = 16)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
