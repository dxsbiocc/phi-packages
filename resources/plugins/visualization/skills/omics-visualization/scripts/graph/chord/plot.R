#!/usr/bin/env Rscript

# Template-ID: graph-chord
#
# Purpose:
#   Draw a weighted Circos chord diagram: sector width is total flow,
#   ribbons are supplied from-to values. This is the non-genomic
#   Circos (migration, trade, cell-cell communication).
#
# Inputs:
#   One row per directed flow. Default example:
#     - from, to: sector names
#     - value: non-negative flow
#
# Output:
#   A PDF, PNG, or SVG chord diagram.
#
# Dependencies:
#   circlize, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Extra sectors are extra names in from/to, not extra ids.
#   directional on/off is CONFIG. Equal-spacing circular networks
#   with node points are graph-circular. Genomic rearrangement
#   ribbons on an ideogram are ideogram-circos. Conserved stacked
#   flows are sankey-*. Do not compute cellchat, BLAST, or a
#   communication score in this script.
#
# Scientific assumptions:
#   from, to, and value are supplied. Sector width is the sum of
#   supplied flows that touch that sector. This script does not
#   symmetricise the matrix unless CONFIG$symmetric is TRUE.

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
        from = "from",
        to = "to",
        value = "value"
    ),
    directional = TRUE,
    symmetric = FALSE,
    palettes = list(
        sector = "Qualitative.Bold"
    ),
    size = list(
        width = 8,
        height = 8
    )
)

load_packages(c("circlize", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
from_col <- config$columns$from
to_col <- config$columns$to
val_col <- config$columns$value
df[[from_col]] <- as.character(df[[from_col]])
df[[to_col]] <- as.character(df[[to_col]])
df[[val_col]] <- as.numeric(df[[val_col]])
if (any(!nzchar(df[[from_col]]) | !nzchar(df[[to_col]]))) {
    stop("from and to must be non-empty.", call. = FALSE)
}
if (any(!is.finite(df[[val_col]]) | df[[val_col]] < 0)) {
    stop("value must be finite and non-negative.", call. = FALSE)
}
if (isTRUE(config$symmetric)) {
    swap <- df
    swap[[from_col]] <- df[[to_col]]
    swap[[to_col]] <- df[[from_col]]
    key <- paste(df[[from_col]], df[[to_col]], sep = "\t")
    rev_key <- paste(swap[[from_col]], swap[[to_col]], sep = "\t")
    df <- rbind(df, swap[!rev_key %in% key, , drop = FALSE])
}

sectors <- unique(c(df[[from_col]], df[[to_col]]))
n_sec <- length(sectors)
sec_cols <- palette_colors(config$palettes$sector, n = n_sec)
names(sec_cols) <- sectors
outline <- palette_colors("Qualitative.Safe")[[12]]

x <- data.frame(
    from = df[[from_col]],
    to = df[[to_col]],
    value = df[[val_col]],
    stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_circos <- function() {
    circlize::circos.par(
        start.degree = 90,
        gap.degree = 3,
        track.margin = c(0.01, 0.01),
        circle.margin = 0.04,
        points.overflow.warning = FALSE
    )
    args <- list(
        x = x,
        grid.col = sec_cols,
        transparency = 0.32,
        annotationTrack = "grid",
        annotationTrackHeight = 0.22,
        link.border = NA,
        grid.border = outline
    )
    if (isTRUE(config$directional)) {
        args$directional <- 1
        args$direction.type <- "diffHeight"
        args$diffHeight <- circlize::mm_h(1.6)
    }
    do.call(circlize::chordDiagram, args)
    circlize::circos.track(
        track.index = 1,
        bg.border = NA,
        panel.fun = function(x, y) {
            sector <- circlize::CELL_META$sector.index
            fill <- unname(sec_cols[[sector]])
            rgb <- grDevices::col2rgb(fill) / 255
            lum <- 0.2126 * rgb[[1]] + 0.7152 * rgb[[2]] + 0.0722 * rgb[[3]]
            txt <- if (lum > 0.45) "black" else "white"
            nch <- nchar(sector)
            cex <- if (nch > 8L) 0.82 else if (nch > 4L) 0.95 else 1.08
            circlize::circos.text(
                circlize::CELL_META$xcenter,
                mean(circlize::CELL_META$ylim),
                sector,
                facing = "bending.inside",
                niceFacing = TRUE,
                adj = c(0.5, 0.5),
                cex = cex,
                font = 2,
                col = txt
            )
        }
    )
    circlize::circos.clear()
}

out_path <- io$output
fmt <- tolower(tools::file_ext(out_path))
if (!fmt %in% c("pdf", "png", "svg")) {
    stop(sprintf("Unsupported output format: %s", fmt), call. = FALSE)
}
out_dir <- dirname(out_path)
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
width <- config$size$width
height <- config$size$height
if (identical(fmt, "png")) {
    grDevices::png(
        out_path, width = width, height = height, units = "in",
        res = 300, bg = "transparent"
    )
} else if (identical(fmt, "pdf")) {
    grDevices::pdf(out_path, width = width, height = height)
} else {
    grDevices::svg(out_path, width = width, height = height, bg = "transparent")
}
graphics::par(mar = c(0.4, 0.4, 0.4, 0.4), bg = "transparent")
tryCatch(
    draw_circos(),
    finally = grDevices::dev.off()
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
