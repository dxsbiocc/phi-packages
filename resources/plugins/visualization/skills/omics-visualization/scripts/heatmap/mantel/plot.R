#!/usr/bin/env Rscript

# Template-ID: heatmap-mantel
#
# Purpose:
#   Draw an upper-triangle correlation heatmap of one variable set,
#   with heart-shaped tiles and curved links to a second variable set.
#
# Inputs:
#   A numeric table (tsv/csv), one row per sample. Default example:
#     - x columns: YAP1, WWTR1
#     - remaining columns: variables shown in the heatmap
#
# Output:
#   A PDF, PNG, or SVG Mantel-style linked correlation heatmap.
#
# Dependencies:
#   ggplot2, readr, linkET
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change which columns are linked or the
#   correlation method and size bins.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Link statistics are pairwise Spearman correlations with cor.test
#   p-values, binned by |r| for line width. Non-significant links
#   (p >= 0.05) are drawn without colour. This is not a Mantel
#   permutation test of distance matrices.

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
        x = c("YAP1", "WWTR1")
    ),
    labels = list(
        title = "",
        x = "",
        y = "",
        fill = "Correlation"
    )
)

load_packages(c("ggplot2", "readr", "linkET"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_cols <- config$columns$x
y_cols <- setdiff(names(df), x_cols)
if (!length(y_cols)) {
    stop("Need at least one column outside columns$x for the heatmap.",
        call. = FALSE
    )
}

for (col in c(x_cols, y_cols)) {
    df[[col]] <- as.numeric(df[[col]])
}

x_mat <- as.data.frame(df[x_cols])
y_mat <- as.data.frame(df[y_cols])

bin_labels <- function(breaks) {
    c(
        paste0("<", breaks[[1]]),
        if (length(breaks) > 1L) {
            paste0(breaks[-length(breaks)], "~", breaks[-1])
        },
        paste0(">=", breaks[[length(breaks)]])
    )
}

r_breaks <- c(0.5, 0.6, 0.7)

mantel_data <- as_md_tbl(suppressMessages(
    correlate(x_mat, y_mat, method = "spearman", cor.test = TRUE)
))
names(mantel_data)[1:2] <- c("x", "y")
r_col <- if ("r" %in% names(mantel_data)) "r" else names(mantel_data)[[3]]
p_col <- if ("p" %in% names(mantel_data)) {
    "p"
} else if ("p.value" %in% names(mantel_data)) {
    "p.value"
} else {
    stop("Correlate table has no p-value column.", call. = FALSE)
}
mantel_data$size <- cut(
    mantel_data[[r_col]],
    breaks = c(-Inf, r_breaks, Inf),
    labels = bin_labels(r_breaks),
    right = FALSE
)
mantel_data$color <- ifelse(
    mantel_data[[p_col]] < 0.05,
    mantel_data[[r_col]],
    NA_real_
)

env_cor <- suppressMessages(
    correlate(y_mat, method = "spearman", cor.test = TRUE)
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_colours <- grDevices::colorRampPalette(c(
    "#2dde98", "#ffffffff", "#ffc168"
))(50)

p <- qcorrplot(env_cor, type = "upper") +
    geom_shaping(marker = marker("heart")) +
    geom_couple(
        data = mantel_data,
        mapping = aes(colour = color, size = size),
        curvature = nice_curvature()
    ) +
    scale_fill_gradientn(
        name = config$labels$fill,
        colours = fill_colours
    ) +
    scale_size_manual(values = c(0.6, 0.9, 1.2, 1.5)) +
    scale_colour_gradientn(
        colours = fill_colours,
        na.value = "transparent"
    ) +
    guides(
        size = guide_legend(
            title = "Spearman's R",
            override.aes = list(colour = "grey50"),
            order = 2
        ),
        fill = guide_colorbar(title = "Spearman's R", order = 3),
        colour = "none"
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme(
        plot.background = element_blank(),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_blank(),
        legend.position = "left",
        panel.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
