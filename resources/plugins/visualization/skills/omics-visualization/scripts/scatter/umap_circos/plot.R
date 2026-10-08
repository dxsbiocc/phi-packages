#!/usr/bin/env Rscript

# Template-ID: scatter-umap-circos
#
# Purpose:
#   Draw a circular composition ring per cell type, with a supplied
#   UMAP embedding in the hole. Outer tracks are stacked proportions
#   of extra categorical columns (time, condition, ...).
#
# Inputs:
#   One row per cell. Default example:
#     - umap_1, umap_2: supplied embedding
#     - Celltype: sector and point colour
#     - Timepoint, Condition: optional composition tracks
#
# Output:
#   A PDF, PNG, or SVG UMAP-in-circos figure.
#
# Dependencies:
#   circlize, ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, tracks, scale, zoom,
#   palette). Discrete hex come from one catalog id; unused colours
#   go to the extra rings so tracks do not repeat celltype hues.
#   Extra metadata rings are extra names in config$tracks, not extra
#   ids. Do not run UMAP, PCA, or clustering. Do not read a Seurat
#   object. Do not wrap plot1cell. A plain embedding is scatter-group.
#   A circos doughnut without embedding is pie-doughnut-bar.
#
# Scientific assumptions:
#   Coordinates and labels are supplied. Sector width follows the
#   number of rows in that cell type (optionally log10 of the cell
#   index). Stacked track segments are within-sector counts, not a
#   statistical test. Point colour matches Celltype.

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
        x = "umap_1",
        y = "umap_2",
        celltype = "Celltype"
    ),
    tracks = c("Timepoint", "Condition"),
    scale = "log10",
    umap_zoom = 0.65,
    point_cex = 0.28,
    point_alpha = 0.35,
    start.degree = 90,
    gap.within = 2,
    gap.last = 2,
    palette = "Qualitative.Prism",
    labels = list(
        celltype = "Celltype"
    ),
    size = list(
        width = 9,
        height = 7
    )
)

load_packages(c("circlize", "ComplexHeatmap", "grid", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
if (length(config$tracks)) {
    require_columns(df, config$tracks)
}

x_col <- config$columns$x
y_col <- config$columns$y
ct_col <- config$columns$celltype
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df[[ct_col]] <- as.character(df[[ct_col]])
if (any(!is.finite(df[[x_col]]) | !is.finite(df[[y_col]]))) {
    stop("umap_1 and umap_2 must be finite.", call. = FALSE)
}
if (any(!nzchar(df[[ct_col]]))) {
    stop("Celltype must be non-empty.", call. = FALSE)
}

ct_n <- sort(table(df[[ct_col]]), decreasing = TRUE)
ct_levels <- names(ct_n)
df[[ct_col]] <- factor(df[[ct_col]], levels = ct_levels)
n_ct <- length(ct_levels)
if (n_ct < 2L) {
    stop("Need at least two cell types for sectors.", call. = FALSE)
}

track_cols <- config$tracks
track_cols <- track_cols[track_cols %in% names(df)]
for (nm in track_cols) {
    df[[nm]] <- factor(as.character(df[[nm]]), levels = unique(as.character(df[[nm]])))
}

df <- df[order(df[[ct_col]]), , drop = FALSE]
df$x_idx <- as.integer(ave(seq_len(nrow(df)), df[[ct_col]], FUN = seq_along))
if (identical(config$scale, "log10")) {
    df$x_polar <- log10(pmax(df$x_idx, 1L))
} else {
    df$x_polar <- as.numeric(df$x_idx)
}

n_track_lv <- vapply(track_cols, function(nm) nlevels(df[[nm]]), integer(1))
pool <- palette_colors(config$palette, n = n_ct + sum(n_track_lv))
ct_cols <- pool[seq_len(n_ct)]
names(ct_cols) <- ct_levels
pool <- pool[-seq_len(n_ct)]
track_cols_map <- list()
for (nm in track_cols) {
    n_lv <- nlevels(df[[nm]])
    cols <- pool[seq_len(n_lv)]
    names(cols) <- levels(df[[nm]])
    track_cols_map[[nm]] <- cols
    pool <- pool[-seq_len(n_lv)]
}

fit_umap <- function(x, y, radius) {
    x <- x - mean(range(x, na.rm = TRUE))
    y <- y - mean(range(y, na.rm = TRUE))
    m <- max(c(abs(x), abs(y), 1e-8))
    list(x = x * radius / m, y = y * radius / m)
}
uv <- fit_umap(df[[x_col]], df[[y_col]], config$umap_zoom)
df$p_x <- uv$x
df$p_y <- uv$y

draw_stack <- function(sector, col_name, fill_map) {
    sub <- df[df[[ct_col]] == sector, , drop = FALSE]
    lv <- levels(df[[col_name]])
    counts <- as.numeric(table(factor(sub[[col_name]], levels = lv)))
    names(counts) <- lv
    counts <- counts[counts > 0]
    if (!length(counts)) {
        return(invisible())
    }
    xlim <- circlize::CELL_META$xlim
    frac <- counts / sum(counts)
    ends <- cumsum(frac)
    starts <- c(0, ends[-length(ends)])
    x0 <- xlim[[1]] + starts * diff(xlim)
    x1 <- xlim[[1]] + ends * diff(xlim)
    circlize::circos.rect(
        x0, 0.15, x1, 0.85,
        col = unname(fill_map[names(counts)]),
        border = NA,
        sector.index = sector
    )
}

axis_ticks <- function(xlim) {
    if (identical(config$scale, "log10")) {
        candidates <- log10(c(1, 10, 100, 1000, 10000))
        at <- candidates[candidates >= xlim[[1]] - 1e-8 &
            candidates <= xlim[[2]] + 1e-8]
        list(
            at = at,
            labels = format(10^at, scientific = FALSE, trim = TRUE)
        )
    } else {
        at <- pretty(xlim, n = 3)
        at <- at[at >= xlim[[1]] - 1e-8 & at <= xlim[[2]] + 1e-8]
        list(
            at = at,
            labels = format(at, scientific = FALSE, trim = TRUE)
        )
    }
}

draw_circos <- function() {
    circlize::circos.clear()
    circlize::circos.par(
        cell.padding = c(0, 0, 0, 0),
        track.margin = c(0.01, 0),
        track.height = 0.01,
        start.degree = config$start.degree,
        gap.degree = c(rep(config$gap.within, n_ct - 1L), config$gap.last),
        canvas.xlim = c(-1.15, 1.55),
        canvas.ylim = c(-1.15, 1.15),
        points.overflow.warning = FALSE
    )
    circlize::circos.initialize(
        sectors = df[[ct_col]],
        x = df$x_polar
    )

    # Outer axis + cell-type names. Same geometry as plot1cell:
    # thin track, x from initialize, circos.axis on the rim, colour as
    # circos.segments at y = 0. Do not put the axis on a ylim=0-1 band.
    circlize::circos.track(
        sectors = df[[ct_col]],
        x = df$x_polar,
        y = rep(0, nrow(df)),
        ylim = c(-1, 1),
        bg.border = NA,
        panel.fun = function(x, y) {
            circlize::circos.text(
                circlize::CELL_META$xcenter,
                circlize::CELL_META$cell.ylim[[2]] + circlize::mm_y(5),
                circlize::CELL_META$sector.index,
                cex = 0.65,
                col = "black",
                facing = "bending.inside",
                niceFacing = TRUE
            )
            ticks <- axis_ticks(circlize::CELL_META$xlim)
            if (length(ticks$at)) {
                circlize::circos.axis(
                    h = "top",
                    major.at = ticks$at,
                    labels = ticks$labels,
                    labels.cex = 0.35,
                    col = "black",
                    labels.col = "black",
                    labels.facing = "inside",
                    labels.niceFacing = TRUE,
                    direction = "outside",
                    minor.ticks = 0,
                    major.tick.length = circlize::mm_y(1)
                )
            }
        }
    )
    for (sector in ct_levels) {
        xr <- range(df$x_polar[df[[ct_col]] == sector])
        circlize::circos.segments(
            x0 = xr[[1]],
            y0 = 0,
            x1 = xr[[2]],
            y1 = 0,
            col = ct_cols[[sector]],
            lwd = 4,
            sector.index = sector
        )
    }

    for (nm in track_cols) {
        fill_map <- track_cols_map[[nm]]
        circlize::circos.track(
            ylim = c(0, 1),
            track.height = 0.045,
            bg.border = NA,
            panel.fun = function(x, y) {
                draw_stack(
                    as.character(circlize::CELL_META$sector.index),
                    nm,
                    fill_map
                )
            }
        )
    }

    graphics::points(
        df$p_x,
        df$p_y,
        pch = 19,
        col = grDevices::adjustcolor(
            unname(ct_cols[as.character(df[[ct_col]])]),
            alpha.f = config$point_alpha
        ),
        cex = config$point_cex
    )
    circlize::circos.clear()

    lgd_list <- list(
        ComplexHeatmap::Legend(
            title = config$labels$celltype,
            at = ct_levels,
            legend_gp = grid::gpar(fill = unname(ct_cols[ct_levels])),
            title_gp = grid::gpar(fontsize = 9),
            labels_gp = grid::gpar(fontsize = 8),
            background = NA
        )
    )
    for (nm in track_cols) {
        lv <- names(track_cols_map[[nm]])
        lgd_list[[length(lgd_list) + 1L]] <- ComplexHeatmap::Legend(
            title = nm,
            at = lv,
            legend_gp = grid::gpar(fill = unname(track_cols_map[[nm]][lv])),
            title_gp = grid::gpar(fontsize = 9),
            labels_gp = grid::gpar(fontsize = 8),
            background = NA
        )
    }
    packed <- ComplexHeatmap::packLegend(list = lgd_list, direction = "vertical")
    ComplexHeatmap::draw(
        packed,
        x = grid::unit(0.92, "npc"),
        y = grid::unit(0.5, "npc"),
        just = "center"
    )
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
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
dpi <- 300
if (identical(fmt, "png")) {
    grDevices::png(
        out_path, width = width, height = height, units = "in",
        res = dpi, bg = "transparent"
    )
} else if (identical(fmt, "pdf")) {
    grDevices::pdf(out_path, width = width, height = height)
} else {
    grDevices::svg(out_path, width = width, height = height, bg = "transparent")
}
old_par <- graphics::par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "transparent")
tryCatch(
    draw_circos(),
    finally = {
        graphics::par(old_par)
        grDevices::dev.off()
    }
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
