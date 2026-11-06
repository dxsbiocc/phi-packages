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
#   directional on/off is CONFIG. A grid of small chords is CONFIG
#   facet (row / col columns): from-sectors take the column colour,
#   other sectors are grey, each unit of value can be its own thin link. Equal-spacing circular networks
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
    # Small-multiple grid: one chord per panel. NULL = single chord.
    # e.g. list(row = c("treatment", "direction"), col = "sender"); the
    # from-sectors of a panel take the column colour, other sectors grey.
    facet = NULL,
    palettes = list(
        sector = "Qualitative.Bold",
        facet_col = "Qualitative.Bold",
        facet_group = "Qualitative.Dark2"
    ),
    receiver_colour = "#d9d9d9",
    facet_label_cex = 0.88,
    facet_unit_links = TRUE,     # draw each unit of value as its own thin link
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
require_columns(df, c(config$columns, as.list(unlist(config$facet))))
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

draw_grid <- function() {
    fr <- config$facet$row
    fc <- config$facet$col
    row_key <- do.call(paste, c(lapply(fr, function(k) as.character(df[[k]])),
        sep = "\r"))
    row_lv <- unique(row_key)
    col_lv <- unique(as.character(df[[fc]]))
    nr <- length(row_lv)
    nc <- length(col_lv)
    col_cols <- setNames(palette_colors(config$palettes$facet_col, n = nc),
        col_lv)
    grp <- vapply(strsplit(row_lv, "\r", fixed = TRUE), `[[`, "", 1)
    sub_lab <- vapply(strsplit(row_lv, "\r", fixed = TRUE),
        function(v) v[[length(v)]], "")
    grp_lv <- unique(grp)
    grp_cols <- setNames(palette_colors(config$palettes$facet_group,
        n = length(grp_lv)), grp_lv)
    two_level <- length(fr) > 1L

    old <- graphics::par(mfrow = c(nr, nc), mar = c(0.2, 0.2, 0.2, 0.2),
        oma = c(0.5, if (two_level) 9 else 6, 2.4, 0.5), xpd = NA)
    on.exit(graphics::par(old), add = TRUE)
    for (i in seq_len(nr)) {
        for (j in seq_len(nc)) {
            sub <- df[row_key == row_lv[[i]] & df[[fc]] == col_lv[[j]], ,
                drop = FALSE]
            if (!nrow(sub)) {
                graphics::plot.new()
                next
            }
            x <- data.frame(from = sub[[from_col]], to = sub[[to_col]],
                value = sub[[val_col]], stringsAsFactors = FALSE)
            if (isTRUE(config$facet_unit_links)) {
                x <- x[rep(seq_len(nrow(x)), pmax(1, round(x$value))), ,
                    drop = FALSE]
                x$value <- 1
            }
            focus <- unique(x$from)
            others <- setdiff(unique(x$to), focus)
            secs <- c(focus, others)
            gc <- setNames(rep(config$receiver_colour, length(secs)), secs)
            gc[focus] <- col_cols[[col_lv[[j]]]]
            circlize::circos.par(start.degree = 90, gap.degree = 3,
                track.margin = c(0.01, 0.01), cell.padding = c(0, 0, 0, 0),
                circle.margin = 0.32, points.overflow.warning = FALSE)
            circlize::chordDiagram(
                x, order = secs, grid.col = gc,
                col = grDevices::adjustcolor(col_cols[[col_lv[[j]]]], 0.55),
                transparency = 0, link.border = "white", link.lwd = 0.4,
                annotationTrack = "grid", annotationTrackHeight = 0.12,
                preAllocateTracks = list(track.height = 0.06)
            )
            circlize::circos.track(track.index = 1, bg.border = NA,
                panel.fun = function(xx, yy) {
                    nm <- circlize::CELL_META$sector.index
                    if (nm %in% focus) return(invisible())
                    circlize::circos.text(circlize::CELL_META$xcenter, 0,
                        nm, facing = "clockwise", niceFacing = TRUE,
                        adj = c(0, 0.5), cex = config$facet_label_cex)
                })
            circlize::circos.clear()
        }
    }
    for (j in seq_len(nc)) {
        graphics::mtext(col_lv[[j]], side = 3, outer = TRUE, line = 0.6,
            at = (j - 0.5) / nc, cex = 0.95)
    }
    for (i in seq_len(nr)) {
        graphics::mtext(sub_lab[[i]], side = 2, outer = TRUE, las = 1,
            line = if (two_level) 0.6 else 0.6, at = 1 - (i - 0.5) / nr,
            adj = 1, cex = 0.85, font = if (two_level) 3 else 1)
    }
    x0 <- graphics::grconvertX(0, "nic", "ndc")
    x1 <- graphics::grconvertX(1, "nic", "ndc")
    y0 <- graphics::grconvertY(0, "nic", "ndc")
    y1 <- graphics::grconvertY(1, "nic", "ndc")
    lab_x <- graphics::grconvertX(0, "nic", "ndc") -
        graphics::grconvertX(0.62 * graphics::par("oma")[[2]], "lines", "ndc") +
        graphics::grconvertX(0, "lines", "ndc")
    graphics::par(mfrow = c(1, 1), fig = c(0, 1, 0, 1), new = TRUE,
        mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
    graphics::rect(x0, y0, x1, y1, border = "grey55", lwd = 0.8)
    if (two_level) {
        runs <- rle(grp)
        ends <- cumsum(runs$lengths)
        starts <- ends - runs$lengths + 1
        for (g in seq_along(runs$values)) {
            yt <- y1 - (starts[[g]] - 1) / nr * (y1 - y0)
            yb <- y1 - ends[[g]] / nr * (y1 - y0)
            if (g > 1) graphics::segments(x0, yt, x1, yt, col = "grey55", lwd = 0.8)
            graphics::text(lab_x, (yt + yb) / 2, runs$values[[g]], font = 2,
                col = grp_cols[[runs$values[[g]]]], adj = c(0.5, 0.5), cex = 1)
        }
    }
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
    if (is.null(config$facet)) draw_circos() else draw_grid(),
    finally = grDevices::dev.off()
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
