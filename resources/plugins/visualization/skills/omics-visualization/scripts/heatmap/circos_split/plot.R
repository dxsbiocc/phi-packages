#!/usr/bin/env Rscript

# Template-ID: heatmap-circos-split
#
# Purpose:
#   Draw a circular heatmap split into sectors by a supplied group,
#   with an optional q-value diamond track and an optional Euler
#   diagram of gene-set overlap in the hole.
#
# Inputs:
#   One row per gene within a group. Default example:
#     - Gene: row label (may repeat across groups)
#     - Group: sector split
#     - qval: supplied adjusted p (diamond fill)
#     - remaining numeric columns: heatmap values (sites / samples)
#
# Output:
#   A PDF, PNG, or SVG split circular heatmap.
#
# Dependencies:
#   circlize, ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, inner, show_qval,
#   start.degree, gaps). inner = "euler" vs "none" and the q-value
#   track are CONFIG, not extra ids. Number of groups or sites is
#   the table, not a new id. A rectangular heatmap with row_split is
#   not a second template. Do not add an id per gene set (Acceleration
#   / Curl / Divergence).
#
# Scientific assumptions:
#   Heatmap cells and qval are supplied. This script does not cluster,
#   rescale, or run a test. Within-group row order is a display sort
#   (qval or file order). The Euler diagram counts gene-name overlap
#   across Group; it is not a statistical test. Euler inset supports
#   2 or 3 groups.

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
        gene = "Gene",
        group = "Group",
        qval = "qval"
    ),
    inner = "euler",
    show_qval = TRUE,
    sort = "qval",
    start.degree = 60,
    gap.within = 5,
    gap.last = 30,
    labels = list(
        fill = "Value",
        qval = "qval",
        inner = "Morphogenic\ngenes"
    ),
    size = list(
        width = 8,
        height = 8
    )
)

load_packages(c("circlize", "ComplexHeatmap", "grid", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

g_gene <- config$columns$gene
g_grp <- config$columns$group
g_q <- config$columns$qval

df[[g_gene]] <- as.character(df[[g_gene]])
df[[g_grp]] <- as.character(df[[g_grp]])
df[[g_q]] <- as.numeric(df[[g_q]])
if (any(!nzchar(df[[g_gene]]) | !nzchar(df[[g_grp]]))) {
    stop("Gene and Group must be non-empty.", call. = FALSE)
}
if (any(!is.finite(df[[g_q]]) | df[[g_q]] < 0 | df[[g_q]] > 1)) {
    stop("qval must be numeric in [0, 1].", call. = FALSE)
}

value_cols <- setdiff(names(df), unlist(config$columns, use.names = FALSE))
if (length(value_cols) < 2L) {
    stop("Need at least two numeric site / sample columns.", call. = FALSE)
}
for (nm in value_cols) {
    df[[nm]] <- as.numeric(df[[nm]])
    if (any(!is.finite(df[[nm]]))) {
        stop(paste("Non-finite value in column", nm), call. = FALSE)
    }
}

grp_levels <- unique(df[[g_grp]])
df[[g_grp]] <- factor(df[[g_grp]], levels = grp_levels)
if (identical(config$sort, "qval")) {
    df <- df[order(df[[g_grp]], df[[g_q]], df[[g_gene]]), , drop = FALSE]
} else {
    df <- df[order(df[[g_grp]]), , drop = FALSE]
}

mat <- as.matrix(df[, value_cols, drop = FALSE])
storage.mode(mat) <- "double"
rownames(mat) <- df[[g_gene]]

heat_pal <- palette_colors("Quantitative.BluGrn")
heat_low <- heat_pal[[1]]
heat_high <- heat_pal[[length(heat_pal)]]
q_pal <- palette_colors("Quantitative.RedOr")
q_light <- q_pal[[1]]
q_dark <- q_pal[[length(q_pal)]]

v_rng <- range(mat, na.rm = TRUE)
if (v_rng[[1]] >= v_rng[[2]]) {
    v_rng[[2]] <- v_rng[[1]] + 1e-6
}
heat_col <- circlize::colorRamp2(v_rng, c(heat_low, heat_high))
# Small qval is darker (more significant).
q_col <- circlize::colorRamp2(c(0, 1), c(q_dark, q_light))

n_grp <- length(grp_levels)
group_cols <- palette_colors("Qualitative.Safe", n = n_grp)
names(group_cols) <- grp_levels

draw_euler <- function(sets, fills, title) {
    n <- length(sets)
    if (n < 2L || n > 3L) {
        grid::grid.text(
            "Euler inset needs 2 or 3 groups",
            gp = grid::gpar(fontsize = 8, col = "#333333")
        )
        return(invisible())
    }
    fills <- unname(fills[names(sets)])
    fills <- grDevices::adjustcolor(fills, alpha.f = 0.35)
    A <- unique(sets[[1]])
    B <- unique(sets[[2]])
    if (n == 2L) {
        circ <- list(
            list(x = 0.38, y = 0.52, r = 0.30),
            list(x = 0.62, y = 0.52, r = 0.30)
        )
        only1 <- setdiff(A, B)
        only2 <- setdiff(B, A)
        both <- intersect(A, B)
        qty <- list(
            list(x = 0.28, y = 0.52, n = length(only1)),
            list(x = 0.72, y = 0.52, n = length(only2)),
            list(x = 0.50, y = 0.52, n = length(both))
        )
    } else {
        C <- unique(sets[[3]])
        circ <- list(
            list(x = 0.38, y = 0.58, r = 0.28),
            list(x = 0.62, y = 0.58, r = 0.28),
            list(x = 0.50, y = 0.36, r = 0.28)
        )
        abc <- Reduce(intersect, list(A, B, C))
        ab <- setdiff(intersect(A, B), C)
        ac <- setdiff(intersect(A, C), B)
        bc <- setdiff(intersect(B, C), A)
        only1 <- setdiff(setdiff(A, B), C)
        only2 <- setdiff(setdiff(B, A), C)
        only3 <- setdiff(setdiff(C, A), B)
        qty <- list(
            list(x = 0.30, y = 0.64, n = length(only1)),
            list(x = 0.70, y = 0.64, n = length(only2)),
            list(x = 0.50, y = 0.26, n = length(only3)),
            list(x = 0.50, y = 0.64, n = length(ab)),
            list(x = 0.36, y = 0.42, n = length(ac)),
            list(x = 0.64, y = 0.42, n = length(bc)),
            list(x = 0.50, y = 0.48, n = length(abc))
        )
    }
    for (i in seq_along(circ)) {
        grid::grid.circle(
            x = circ[[i]]$x,
            y = circ[[i]]$y,
            r = circ[[i]]$r,
            default.units = "npc",
            gp = grid::gpar(
                fill = fills[[i]],
                col = "black",
                lty = 2,
                lwd = 1.2
            )
        )
    }
    for (q in qty) {
        if (q$n <= 0) {
            next
        }
        grid::grid.text(
            as.character(q$n),
            x = q$x,
            y = q$y,
            default.units = "npc",
            gp = grid::gpar(fontsize = 8, col = "#222222")
        )
    }
    if (nzchar(title)) {
        grid::grid.text(
            title,
            x = 0.5,
            y = -0.04,
            default.units = "npc",
            gp = grid::gpar(fontsize = 9, lineheight = 0.9, col = "#222222")
        )
    }
}

draw_circos <- function() {
    circlize::circos.par(
        start.degree = config$start.degree,
        gap.after = c(rep(config$gap.within, n_grp - 1L), config$gap.last),
        track.margin = c(0, 0.01),
        cell.padding = c(0, 0, 0, 0),
        circle.margin = 0.02
    )
    circlize::circos.heatmap(
        mat,
        split = df[[g_grp]],
        cluster = FALSE,
        bg.border = "black",
        bg.lwd = 1,
        cell.border = "white",
        cell.lwd = 0.4,
        rownames.side = "outside",
        rownames.cex = 0.55,
        col = heat_col,
        track.height = 0.28
    )
    heat_track <- circlize::get.current.track.index()
    circlize::circos.track(
        track.index = heat_track,
        bg.border = NA,
        panel.fun = function(x, y) {
            if (circlize::CELL_META$sector.numeric.index != n_grp) {
                return(invisible())
            }
            cn <- colnames(mat)
            n <- length(cn)
            cell_h <- diff(circlize::CELL_META$cell.ylim) / n
            y_coords <- seq(
                circlize::CELL_META$cell.ylim[[1]] + cell_h / 2,
                circlize::CELL_META$cell.ylim[[2]] - cell_h / 2,
                length.out = n
            )
            x0 <- circlize::CELL_META$cell.xlim[[2]]
            x1 <- x0 + circlize::convert_x(1, "mm")
            for (i in seq_len(n)) {
                circlize::circos.lines(
                    c(x0, x1),
                    c(y_coords[[i]], y_coords[[i]]),
                    col = "black",
                    lwd = 1
                )
            }
            circlize::circos.text(
                rep(x0, n) + circlize::convert_x(1.6, "mm"),
                y_coords,
                cn,
                cex = 0.55,
                adj = c(0, 0.5),
                facing = "inside",
                col = "black"
            )
        }
    )

    if (isTRUE(config$show_qval)) {
        circlize::circos.track(
            ylim = c(0, 1),
            track.height = 0.06,
            bg.border = NA,
            panel.fun = function(x, y) {
                sector <- as.character(circlize::CELL_META$sector.index)
                sub <- df[df[[g_grp]] == sector, , drop = FALSE]
                n <- nrow(sub)
                if (!n) {
                    return(invisible())
                }
                x_at <- circlize::CELL_META$xlim[[1]] +
                    (circlize::CELL_META$xlim[[2]] - circlize::CELL_META$xlim[[1]]) *
                    (seq_len(n) - 0.5) / n
                circlize::circos.points(
                    x_at,
                    rep(0.5, n),
                    pch = 18,
                    cex = 0.9,
                    col = q_col(sub[[g_q]])
                )
                if (circlize::CELL_META$sector.numeric.index == n_grp) {
                    x0 <- circlize::CELL_META$cell.xlim[[2]]
                    circlize::circos.lines(
                        c(x0, x0 + circlize::convert_x(1, "mm")),
                        c(0.5, 0.5),
                        col = "black",
                        lwd = 1
                    )
                    circlize::circos.text(
                        x0 + circlize::convert_x(1.6, "mm"),
                        0.5,
                        config$labels$qval,
                        cex = 0.55,
                        adj = c(0, 0.5),
                        facing = "inside",
                        col = "black"
                    )
                }
            }
        )
    }

    circlize::circos.track(
        ylim = c(0, 1),
        track.height = 0.07,
        bg.col = grDevices::adjustcolor(group_cols[grp_levels], alpha.f = 0.35),
        bg.border = NA,
        panel.fun = function(x, y) {
            circlize::circos.text(
                circlize::CELL_META$xcenter,
                0.45,
                circlize::CELL_META$sector.index,
                facing = "bending.inside",
                niceFacing = TRUE,
                cex = 0.7,
                adj = c(0.5, 0.5),
                col = "black"
            )
        }
    )
    circlize::circos.clear()

    heat_lgd <- ComplexHeatmap::Legend(
        title = config$labels$fill,
        col_fun = heat_col,
        at = pretty(v_rng, n = 4),
        title_position = "leftcenter-rot",
        title_gp = grid::gpar(fontsize = 9),
        labels_gp = grid::gpar(fontsize = 8),
        background = NA
    )
    ComplexHeatmap::draw(
        heat_lgd,
        x = grid::unit(0.98, "npc") - grid::unit(2, "mm"),
        y = grid::unit(0.92, "npc"),
        just = c("right", "top")
    )
    if (isTRUE(config$show_qval)) {
        q_lgd <- ComplexHeatmap::Legend(
            title = config$labels$qval,
            col_fun = q_col,
            at = c(0, 0.25, 0.5, 0.75, 1),
            title_position = "leftcenter-rot",
            title_gp = grid::gpar(fontsize = 9),
            labels_gp = grid::gpar(fontsize = 8),
            background = NA
        )
        ComplexHeatmap::draw(
            q_lgd,
            x = grid::unit(0.98, "npc") - grid::unit(2, "mm"),
            y = grid::unit(0.22, "npc"),
            just = c("right", "bottom")
        )
    }

    if (identical(config$inner, "euler")) {
        sets <- lapply(grp_levels, function(g) {
            unique(df[[g_gene]][df[[g_grp]] == g])
        })
        names(sets) <- grp_levels
        grid::pushViewport(grid::viewport(
            x = 0.5,
            y = 0.5,
            width = 0.26,
            height = 0.26
        ))
        draw_euler(sets, group_cols, config$labels$inner)
        grid::popViewport()
    }
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
