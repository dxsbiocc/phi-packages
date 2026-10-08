#!/usr/bin/env Rscript

# Template-ID: heatmap-mantel-size
#
# Purpose:
#   Draw an upper-triangle correlation heatmap of one variable set, with
#   sized and coloured links to a second (Mantel-style) variable set.
#
# Inputs:
#   A numeric table (tsv/csv), one row per sample. Default example:
#     - x columns: YAP1, WWTR1
#     - remaining columns: environment / omics variables in the heatmap
#
# Output:
#   A PDF, PNG, or SVG Mantel-style linked correlation heatmap.
#
# Dependencies:
#   ggplot2, readr, ggcor, linkET
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change which columns are linked or the
#   correlation method and size/p-value bins.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Link statistics are pairwise Pearson correlations with cor.test
#   p-values, displayed in Mantel-style bins; this is not a Mantel
#   permutation test of distance matrices.
#   Heatmap tiles are Pearson correlations among the non-x columns.

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

load_packages(c("ggplot2", "readr", "ggcor"))

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
p_breaks <- c(0.001, 0.01, 0.05)

mantel_data <- linkET::as_md_tbl(suppressWarnings(
    correlate(x_mat, y_mat, method = "pearson", cor.test = TRUE)
))
names(mantel_data)[1:2] <- c("spec", "env")
r_col <- if ("r" %in% names(mantel_data)) "r" else names(mantel_data)[[3]]
p_col <- if ("p.value" %in% names(mantel_data)) {
    "p.value"
} else if ("p" %in% names(mantel_data)) {
    "p"
} else {
    stop("Correlate table has no p-value column.", call. = FALSE)
}
mantel_data$size <- cut(
    mantel_data[[r_col]],
    breaks = c(-Inf, r_breaks, Inf),
    labels = bin_labels(r_breaks),
    right = FALSE
)
mantel_data$p.value <- cut(
    mantel_data[[p_col]],
    breaks = c(-Inf, p_breaks, Inf),
    labels = bin_labels(p_breaks),
    right = FALSE
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_colours <- grDevices::colorRampPalette(c(
    "#87c55f", "white", "#fe88b1"
))(50)

p <- quickcor(y_mat, type = "upper") +
    geom_square() +
    add_link(
        mantel_data,
        mapping = aes(colour = p.value, size = size)
    ) +
    scale_fill_gradientn(
        name = config$labels$fill,
        colours = fill_colours
    ) +
    scale_size_manual(values = c(0.2, 0.4, 0.7, 1.1)) +
    scale_colour_manual(values = c(
        "#1aafd0", "#3be8b0", "#6a67ce", "#caccd1"
    )) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme(
        plot.background = element_blank(),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_blank(),
        legend.position = "right",
        panel.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
