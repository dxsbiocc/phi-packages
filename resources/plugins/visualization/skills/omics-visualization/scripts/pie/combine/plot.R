#!/usr/bin/env Rscript

# Template-ID: pie-combine
#
# Purpose:
#   Draw a two-ring pie: inner ring for group totals, outer ring for items.
#
# Inputs:
#   A table with one row per item. Default example:
#     - name: item label
#     - value: numeric size
#     - category: group used for the inner ring
#
# Output:
#   A PDF, PNG, or SVG combined nested pie chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or reshape the table.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Inner-ring values are sums of items within each category.
#   Outer-ring slices are the original item values.
#   Values are treated as already summarized, not as raw replicates.

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
        group = "category"
    ),
    labels = list(
        title = "Combined Pie Chart",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group
df[[y_col]] <- as.numeric(df[[y_col]])
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
inner <- stats::aggregate(
    df[[y_col]],
    by = list(group = df[[group_col]]),
    FUN = sum
)
names(inner)[2] <- "value"
inner$ring_x <- 1
inner$fill <- as.character(inner$group)
outer <- df
outer$ring_x <- 2
outer$fill <- as.character(outer[[x_col]])
plot_df <- rbind(
    data.frame(
        ring_x = inner$ring_x,
        value = inner$value,
        fill = inner$fill,
        stringsAsFactors = FALSE
    ),
    data.frame(
        ring_x = outer$ring_x,
        value = outer[[y_col]],
        fill = outer$fill,
        stringsAsFactors = FALSE
    )
)
plot_df$fill <- factor(plot_df$fill, levels = unique(plot_df$fill))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(plot_df, aes(x = ring_x, y = value, fill = fill)) +
    geom_col(width = 1, colour = "white") +
    scale_fill_manual(
        values = expand_palette(
            c("#0079bf",
            "#70b500",
            "#ff9f1a",
            "#eb5a46",
            "#f2d600",
            "#c377e0",
            "#ff78cb",
            "#00c2e0",
            "#51e898"),
            nlevels(plot_df$fill)
        )
    ) +
    coord_polar(theta = "y") +
    xlim(0, 2.6) +
    labs(title = config$labels$title, x = NULL, y = NULL, fill = "") +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 8)
