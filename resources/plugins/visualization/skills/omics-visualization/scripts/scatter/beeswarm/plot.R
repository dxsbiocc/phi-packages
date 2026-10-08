#!/usr/bin/env Rscript

# Template-ID: scatter-beeswarm
#
# Purpose:
#   Draw a beeswarm plot of a numeric value across categories, with
#   quartile bars and median markers per group.
#
# Inputs:
#   A table with one row per observation. Default example:
#     - name: category label
#     - value: numeric value
#     - group: grouping label
#
# Output:
#   A PDF, PNG, or SVG beeswarm plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggbeeswarm
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or the quartile
#   summary.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one observation, not a pre-summarized statistic.
#   Quartile bars and the median marker are descriptive, not a test.
#   Category order is a display choice.

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
        x = "name",
        y = "value",
        group = "group"
    ),
    labels = list(
        title = "",
        x = "",
        y = "Expression level (log2 TPM)"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggbeeswarm"))

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

summary_data <- do.call(rbind, lapply(
    split(df, list(df[[x_col]], df[[group_col]]), drop = TRUE),
    function(sub) {
        v <- sub[[y_col]]
        data.frame(
            name = sub[[x_col]][[1]],
            group = sub[[group_col]][[1]],
            first = as.numeric(stats::quantile(v, 0.25)),
            second = as.numeric(stats::quantile(v, 0.5)),
            third = as.numeric(stats::quantile(v, 0.75))
        )
    }
))
names(summary_data)[1:2] <- c(x_col, group_col)
summary_data[[x_col]] <- factor(summary_data[[x_col]], levels = levels(df[[x_col]]))
summary_data[[group_col]] <- factor(
    summary_data[[group_col]],
    levels = levels(df[[group_col]])
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c("#f784b6", "#79ceb8")

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]],
    colour = .data[[config$columns$group]]
)) +
    geom_quasirandom(width = 0.15, size = 3, dodge.width = 0.8, alpha = 0.6) +
    geom_errorbar(
        data = summary_data,
        aes(
            x = .data[[config$columns$x]],
            ymin = first,
            ymax = third,
            group = .data[[config$columns$group]]
        ),
        width = 0,
        position = position_dodge(width = 0.8),
        linewidth = 1,
        inherit.aes = FALSE
    ) +
    geom_point(
        data = summary_data,
        aes(
            x = .data[[config$columns$x]],
            y = second,
            group = .data[[config$columns$group]]
        ),
        position = position_dodge2(width = 0.8),
        size = 5,
        shape = 21,
        fill = "white",
        stroke = 1,
        inherit.aes = FALSE
    ) +
    geom_point(
        data = summary_data,
        aes(
            x = .data[[config$columns$x]],
            y = second,
            fill = .data[[config$columns$group]],
            group = .data[[config$columns$group]]
        ),
        position = position_dodge2(width = 0.8),
        size = 4,
        shape = 21,
        stroke = 0,
        inherit.aes = FALSE
    ) +
    scale_colour_manual(values = colour_values) +
    scale_fill_manual(values = colour_values) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.8),
        legend.background = element_rect(fill = "white", colour = "grey80"),
        axis.line = element_blank(),
        panel.border = element_rect(fill = "transparent", colour = "black"),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
