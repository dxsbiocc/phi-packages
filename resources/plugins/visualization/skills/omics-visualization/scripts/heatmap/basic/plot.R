#!/usr/bin/env Rscript

# Template-ID: heatmap-basic
#
# Purpose:
#   Draw a tile heatmap of a numeric value across hours of day and days of week.
#
# Inputs:
#   A table with one id column and numeric day columns. Default example:
#     - hours: row label (time of day)
#     - remaining columns: one numeric series per weekday
#
# Output:
#   A PDF, PNG, or SVG heatmap.
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
#   Each cell is one already summarized count or intensity.
#   Row and column order in the file are display order, not a clustering result.
#   No scaling or clustering is applied.

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
        id = "hours"
    ),
    labels = list(
        title = "Hourly Heatmap",
        x = "Day",
        y = "Hour",
        fill = "Count"
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
id_col <- config$columns$id
value_cols <- setdiff(names(df), id_col)
if (!length(value_cols)) {
    stop("Heatmap input needs at least one numeric column besides the id.",
        call. = FALSE
    )
}
for (col in value_cols) {
    df[[col]] <- as.numeric(df[[col]])
}
df[[id_col]] <- factor(df[[id_col]], levels = unique(df[[id_col]]))
long <- data.frame(
    id = rep(df[[id_col]], times = length(value_cols)),
    column = factor(
        rep(value_cols, each = nrow(df)),
        levels = value_cols
    ),
    value = as.numeric(as.matrix(df[value_cols])),
    stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(long, aes(x = column, y = id, fill = value)) +
    geom_tile(colour = "white") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(levels(long$id)), expand = c(0, 0)) +
    scale_fill_gradientn(colours = c("#cfd2d7", "#a8dab0", "#1d4f60")) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        fill = config$labels$fill
    ) +
    theme_prism() +
    theme(
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
