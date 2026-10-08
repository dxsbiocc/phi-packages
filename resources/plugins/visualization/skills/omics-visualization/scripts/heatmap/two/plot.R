#!/usr/bin/env Rscript

# Template-ID: heatmap-two
#
# Purpose:
#   Draw a two-set correlation heatmap comparing one group of variables
#   against another (stars for correlation, crosses for non-significance).
#
# Inputs:
#   A numeric table (tsv/csv), one row per sample and one column per
#   variable. Default example:
#     - x columns: YAP1, WWTR1, TEAD1–4, SMAD1–4
#     - remaining columns: the second variable set
#
# Output:
#   A PDF, PNG, or SVG two-set correlation heatmap.
#
# Dependencies:
#   ggplot2, readr, ggcor
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change which columns form each set or
#   the correlation method.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Rows are paired observations used to compute Spearman correlations.
#   Crosses mark pairs with p >= 0.05 from cor.test; they are a display
#   overlay, not a multiple-testing correction.

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
        x = c(
            "YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4",
            "SMAD1", "SMAD2", "SMAD3", "SMAD4"
        )
    ),
    labels = list(
        title = "",
        x = "",
        y = "",
        fill = "Correlation"
    )
)

load_packages(c("ggplot2", "readr", "ggcor"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_cols <- config$columns$x
y_cols <- setdiff(names(df), x_cols)
if (!length(y_cols)) {
    stop("Need at least one column outside columns$x for the second set.",
        call. = FALSE
    )
}

for (col in c(x_cols, y_cols)) {
    df[[col]] <- as.numeric(df[[col]])
}

cor_tbl <- as_cor_tbl(correlate(
    as.data.frame(df[x_cols]),
    as.data.frame(df[y_cols]),
    method = "spearman",
    cor.test = TRUE
))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_colours <- grDevices::colorRampPalette(c(
    "#009dd3", "#ffffffff", "#f29c98"
))(50)

p <- quickcor(cor_tbl) +
    geom_star() +
    geom_cross(colour = "#f5e797") +
    scale_fill_gradientn(
        name = config$labels$fill,
        colours = fill_colours
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme(
        plot.background = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        panel.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
