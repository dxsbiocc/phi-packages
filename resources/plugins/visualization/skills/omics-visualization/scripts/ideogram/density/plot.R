#!/usr/bin/env Rscript

# Template-ID: ideogram-density
#
# Purpose:
#   Draw chromosomes as ideogram tracks filled by a supplied continuous
#   value (gene density, FST, etc.).
#
# Inputs:
#   One row per window. Default example:
#     - Chrom: chromosome name
#     - Start, End: window coordinates
#     - Value: continuous statistic
#
# Output:
#   A PDF, PNG, or SVG density ideogram.
#
# Dependencies:
#   ggplot2, readr, ggideogram (GitHub: dxsbiocc/ggideogram)
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, flip, gradient).
#   G-banding is ideogram-karyotype. Transcript exon models are
#   ideogram-gene. A second density track (LTR vs genes) is a second
#   geom_ideogram layer in PLOT, not a new id. RNA points and gene
#   labels are overlays, not extra templates.
#
# Scientific assumptions:
#   Value is a supplied window statistic, not computed here.
#   Windows are already binned. Fill is a two-end catalog gradient.

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
        chrom = "Chrom",
        start = "Start",
        end = "End",
        value = "Value"
    ),
    flip = FALSE,
    radius_pt = 4,
    width = 0.5,
    labels = list(
        title = "",
        fill = "Density"
    )
)

load_packages(c("ggplot2", "readr", "ggideogram"))
patch_ggideogram_ggplot2()

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

ch <- config$columns$chrom
st <- config$columns$start
en <- config$columns$end
vl <- config$columns$value

df[[ch]] <- factor(df[[ch]], levels = unique(df[[ch]]))
df[[st]] <- as.numeric(df[[st]])
df[[en]] <- as.numeric(df[[en]])
df[[vl]] <- as.numeric(df[[vl]])

grad <- palette_colors("Quantitative.BluGrn")
fill_low <- grad[[1]]
fill_high <- grad[[length(grad)]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot2::ggplot(df) +
    ggideogram::geom_ideogram(
        ggplot2::aes(
            x = .data[[ch]],
            ymin = .data[[st]],
            ymax = .data[[en]],
            chrom = .data[[ch]],
            fill = .data[[vl]]
        ),
        width = config$width,
        radius = grid::unit(config$radius_pt, "pt"),
        linewidth = 0.4,
        chrom.col = "#888888",
        show.legend = TRUE
    ) +
    ggplot2::scale_fill_gradient(
        low = fill_low,
        high = fill_high,
        name = config$labels$fill
    ) +
    ggplot2::labs(title = config$labels$title, x = NULL, y = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(
        plot.background = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(colour = "black", size = 9),
        plot.margin = ggplot2::margin(8, 8, 8, 8)
    )

if (isTRUE(config$flip)) {
    p <- p +
        ggplot2::coord_flip() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(colour = "black", size = 9)
        )
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 5.5)
