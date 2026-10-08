#!/usr/bin/env Rscript

# Template-ID: ideogram-karyotype
#
# Purpose:
#   Draw chromosome ideograms with G-banding (stain) and p/q arms.
#
# Inputs:
#   One row per cytoband. Default example:
#     - Chrom: chromosome name
#     - Start, End: band coordinates
#     - Stain: gneg / gpos* / acen / gvar / stalk
#     - Arm: p or q
#
# Output:
#   A PDF, PNG, or SVG karyotype ideogram.
#
# Dependencies:
#   ggplot2, readr, ggideogram (GitHub: dxsbiocc/ggideogram)
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, flip, which
#   chromosomes). coord_flip and hg19 vs mm10 are CONFIG / data, not
#   extra ids. Gene models are ideogram-gene. Continuous density fill
#   is ideogram-density. A supplied gene-loci table with labels is
#   ideogram-loci (circular vs vertical is CONFIG there).
#
# Scientific assumptions:
#   Bands and stains are supplied (UCSC cytoBand or equivalent). This
#   script does not fetch a genome. Fill uses ggideogram cytoband
#   colours (standard G-band encoding).

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
        stain = "Stain",
        arm = "Arm"
    ),
    flip = FALSE,
    radius_pt = 5,
    width = 0.45,
    labels = list(
        title = ""
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
sn <- config$columns$stain
ar <- config$columns$arm

df[[ch]] <- factor(df[[ch]], levels = unique(df[[ch]]))
df[[st]] <- as.numeric(df[[st]])
df[[en]] <- as.numeric(df[[en]])
df[[sn]] <- as.character(df[[sn]])
df[[ar]] <- as.character(df[[ar]])

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
            fill = .data[[sn]],
            arm = .data[[ar]]
        ),
        width = config$width,
        radius = grid::unit(config$radius_pt, "pt"),
        linewidth = 0.2,
        colour = "black",
        show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
        values = ggideogram::cytoband_colors,
        na.value = "white"
    ) +
    ggplot2::labs(title = config$labels$title, x = NULL, y = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(
        plot.background = ggplot2::element_blank(),
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
save_ggplot(p, io$output, width = 8, height = 6)
