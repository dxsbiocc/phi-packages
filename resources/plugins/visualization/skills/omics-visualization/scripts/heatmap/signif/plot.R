#!/usr/bin/env Rscript

# Template-ID: heatmap-signif
#
# Purpose:
#   Draw a correlation heatmap with significance overlays: squares and
#   shade on the lower triangle, circles and marks on the upper triangle,
#   and crosses on weak correlations.
#
# Inputs:
#   A numeric table (tsv/csv), one row per sample and one column per
#   variable. Default example uses Hippo-pathway genes.
#
# Output:
#   A PDF, PNG, or SVG correlation heatmap with significance marks.
#
# Dependencies:
#   ggplot2, readr, ggcor
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change the correlation method or which
#   columns are included.
#   Edit PLOT only when overlay filters or styling must change.
#
# Scientific assumptions:
#   Pairwise Spearman correlations and cor.test p-values are computed
#   from paired rows. Marks and crosses are display filters
#   (p < 0.05, |r| thresholds), not a multiple-testing correction.

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
    columns = list(),
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

feature_cols <- names(df)
if (length(feature_cols) < 2L) {
    stop("Need at least two numeric columns to compute correlations.",
        call. = FALSE
    )
}
for (col in feature_cols) {
    df[[col]] <- as.numeric(df[[col]])
}

cor_tbl <- as_cor_tbl(correlate(
    as.data.frame(df[feature_cols]),
    method = "spearman",
    cor.test = TRUE
))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_colours <- grDevices::colorRampPalette(c(
    "#1cc7d0", "#ffffffff", "#fbb05b"
))(50)

p <- quickcor(cor_tbl) +
    geom_square(
        mapping = aes(fill = r),
        data = get_data(type = "lower")
    ) +
    geom_shade(data = get_data(type = "lower"), sign = 1) +
    geom_circle2(
        mapping = aes(fill = r),
        data = get_data(type = "upper")
    ) +
    geom_mark(
        data = get_data(
            r != 1, p.value < 0.05, abs(r) >= 0.3,
            type = "upper"
        ),
        size = 3,
        sep = "\n"
    ) +
    geom_cross(data = get_data(abs(r) < 0.5)) +
    geom_diag_label(size = 2.5) +
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
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        panel.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
