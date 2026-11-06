#!/usr/bin/env Rscript

# Template-ID: heatmap-circos-cohort
#
# Purpose:
#   Circular per-sample cohort overview: one equal slice per sample
#   with its ID outside, numeric lollipop tracks (with a y axis) and
#   categorical tile tracks, plus summary columns in the gap aligned
#   with each ring (box + dots with quartile labels for numeric
#   tracks, percent stacked bar for categorical tracks) and legends.
#   Summary positions come from the ring radii, so nothing is placed
#   by hand.
#
# Inputs:
#   One row per sample. Default example:
#     - ID: sample id (slice order = file order unless sort_by set)
#     - Age, BMI: numeric tracks (NA = no lollipop)
#     - Gender, Ethnicity, IRIS: categorical tracks (NA = "NA")
#
# Output:
#   A PDF, PNG, or SVG circular cohort figure.
#
# Dependencies:
#   circlize, ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG: one `tracks` entry per column, outer to inner
#   (type lollipop / tile, colour or per-level colours, height,
#   summary title). Extra attributes are extra entries; gap, start
#   angle, ID labels, axis and quartile labels are CONFIG. Grouped
#   category sectors with totals are heatmap-circos-annotation.
#
# Scientific assumptions:
#   Values are supplied. Quartiles, medians and percentages in the
#   gap are descriptive display summaries of the supplied values;
#   no test is run. Percentages use non-missing samples.

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
    id = "ID",
    sort_by = NULL,              # NULL = file order; or a column name
    # outer to inner; colours are "Family.Name:index"
    tracks = list(
        list(column = "Age", type = "lollipop", color = "Qualitative.Bold:1",
            height = 0.15, title = "Age\n(years)"),
        list(column = "BMI", type = "lollipop", color = "Qualitative.Vivid / 10:10",
            height = 0.15, title = "BMI"),
        list(column = "Gender", type = "tile", height = 0.065,
            colors = c(Female = "Qualitative.Bold:5",
                Male = "Qualitative.Bold:8"),
            title = "Gender (%)"),
        list(column = "Ethnicity", type = "tile", height = 0.065,
            palette = "Qualitative.Pastel", title = "Ethnicity (%)"),
        list(column = "IRIS", type = "tile", height = 0.065,
            colors = c(IR = "Qualitative.Vivid:2", IS = "Qualitative.Vivid / 10:10",
                Unknown = "Qualitative.Safe:12"),
            title = "Insulin\nresponse (%)", legend = "Insulin response")
    ),
    gap = 90,                    # degrees left empty for the summaries
    start_degree = 90,
    show_ids = TRUE,
    id_track_height = 0.08,
    id_cex = 0.42,
    show_axis = TRUE,
    quartile_labels = TRUE,
    summary_height = 0.62,       # summary column height, canvas units
    summary_offset = 0.04,       # gap between ring end and column
    legend_xy = c(0.02, 0.98),   # npc, top-left anchor
    size = list(width = 10, height = 10)
)

load_packages(c("circlize", "ComplexHeatmap", "grid", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
cols_needed <- c(config$id, vapply(config$tracks, `[[`, "", "column"))
require_columns(df, as.list(stats::setNames(cols_needed, cols_needed)))
df[[config$id]] <- as.character(df[[config$id]])
if (anyDuplicated(df[[config$id]])) {
    stop("Sample IDs must be unique.", call. = FALSE)
}
if (!is.null(config$sort_by)) {
    df <- df[order(df[[config$sort_by]]), , drop = FALSE]
}
ids <- df[[config$id]]
n <- length(ids)

trk <- lapply(config$tracks, function(t) {
    v <- df[[t$column]]
    if (identical(t$type, "lollipop")) {
        t$v <- as.numeric(v)
        t$col <- pref(t$color)
    } else {
        v <- as.character(v)
        v[is.na(v) | !nzchar(v)] <- "NA"
        lv <- unique(v)
        lv <- c(setdiff(lv, c("Unknown", "NA")), intersect(c("Unknown", "NA"),
            lv))
        if (!is.null(t$colors)) {
            lv <- c(names(t$colors)[names(t$colors) %in% lv],
                setdiff(lv, names(t$colors)))
        }
        fall <- palette_colors(if (is.null(t$palette)) "Qualitative.Safe"
            else t$palette, n = length(lv))
        t$pal <- stats::setNames(vapply(seq_along(lv), function(i) {
            if (!is.null(t$colors) && lv[[i]] %in% names(t$colors)) {
                pref(t$colors[[lv[[i]]]])
            } else {
                fall[[i]]
            }
        }, character(1)), lv)
        t$v <- v
        t$lv <- lv
    }
    t
})

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_figure <- function() {
    circlize::circos.par(
        start.degree = config$start_degree, clock.wise = TRUE,
        gap.after = c(rep(0, n - 1), config$gap),
        cell.padding = c(0, 0, 0, 0), track.margin = c(0.004, 0.004),
        circle.margin = c(0.02, 0.02, 0.02, 0.02),
        points.overflow.warning = FALSE
    )
    circlize::circos.initialize(sectors = factor(ids, levels = ids),
        xlim = c(0.5, 1.5))
    if (isTRUE(config$show_ids)) {
        circlize::circos.track(ylim = c(0, 1), bg.border = NA,
            track.height = config$id_track_height,
            panel.fun = function(x, y) {
                circlize::circos.text(1, 0, circlize::CELL_META$sector.index,
                    facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.5), cex = config$id_cex)
            })
    }
    first <- ids[[1]]
    last <- ids[[n]]
    band <- list()
    for (k in seq_along(trk)) {
        t <- trk[[k]]
        if (identical(t$type, "lollipop")) {
            rng <- range(t$v, na.rm = TRUE)
            ylim <- c(rng[[1]] - diff(rng) * 0.08, rng[[2]] + diff(rng) * 0.08)
            circlize::circos.track(ylim = ylim, bg.border = "black",
                bg.lwd = 0.5, track.height = t$height,
                panel.fun = function(x, y) {
                    i <- circlize::CELL_META$sector.numeric.index
                    v <- t$v[[i]]
                    if (is.finite(v)) {
                        circlize::circos.lines(c(1, 1), c(ylim[[1]], v),
                            col = t$col, lwd = 1.4)
                        circlize::circos.points(1, v, pch = 16, cex = 0.55,
                            col = t$col)
                    }
                    if (isTRUE(config$show_axis) && i == 1L) {
                        at <- round(c(rng[[1]], mean(rng), rng[[2]]))
                        circlize::circos.yaxis(side = "left", at = at,
                            labels.cex = 0.55, labels.niceFacing = FALSE,
                            tick.length = circlize::convert_x(0.6, "mm"))
                    }
                })
        } else {
            circlize::circos.track(ylim = c(0, 1), bg.border = "black",
                bg.lwd = 0.5, track.height = t$height,
                panel.fun = function(x, y) {
                    i <- circlize::CELL_META$sector.numeric.index
                    circlize::circos.rect(0.5, 0, 1.5, 1,
                        col = t$pal[[t$v[[i]]]], border = "black", lwd = 0.3)
                })
        }
        band[[k]] <- c(
            circlize::get.cell.meta.data("cell.bottom.radius", last,
                circlize::get.current.track.index()),
            circlize::get.cell.meta.data("cell.top.radius", last,
                circlize::get.current.track.index()))
    }
    end_deg <- circlize::get.cell.meta.data("cell.end.degree", last, 1)
    circlize::circos.clear()

    # Summary columns: above the ring ends at the gap boundary.
    graphics::par(xpd = NA)
    ang <- end_deg / 180 * pi
    ux <- c(cos(ang), sin(ang))          # radial direction at the boundary
    nx <- c(-ux[[2]], ux[[1]])           # perpendicular, into the gap
    if (sum(nx * c(cos(ang - 0.1), sin(ang - 0.1))) < 0) nx <- -nx
    txt_srt <- atan2(nx[[2]], nx[[1]]) * 180 / pi - 90
    H <- config$summary_height
    off <- config$summary_offset
    to_xy <- function(r, h) cbind(r * ux[[1]] + h * nx[[1]],
        r * ux[[2]] + h * nx[[2]])
    quad <- function(r0, r1, h0, h1, ...) {
        p <- rbind(to_xy(r0, h0), to_xy(r1, h0), to_xy(r1, h1), to_xy(r0, h1))
        graphics::polygon(p[, 1], p[, 2], ...)
    }
    for (k in seq_along(trk)) {
        t <- trk[[k]]
        r0 <- band[[k]][[1]]
        r1 <- band[[k]][[2]]
        w <- r1 - r0
        quad(r0, r1, off, off + H, border = "black", lwd = 0.6)
        if (identical(t$type, "lollipop")) {
            v <- t$v[is.finite(t$v)]
            rng <- range(v)
            hy <- function(z) off + H * (0.06 + 0.88 * (z - rng[[1]]) /
                max(diff(rng), 1e-9))
            q <- stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
            bins <- cut(v, 30, labels = FALSE)
            rk <- stats::ave(seq_along(v), bins, FUN = seq_along)
            sgn <- ifelse(rk %% 2 == 0, 1, -1)
            jit <- sgn * ceiling((rk - 1) / 2) * w * 0.07
            jit <- pmax(pmin(jit, w * 0.42), -w * 0.42)
            rc <- (r0 + r1) / 2
            xy <- to_xy(rc + jit, hy(v))
            graphics::points(xy[, 1], xy[, 2], pch = 16, cex = 0.45,
                col = grDevices::adjustcolor(t$col, 0.8))
            seg <- function(ra, rb, h, lwd = 1) {
                p <- rbind(to_xy(ra, h), to_xy(rb, h))
                graphics::lines(p[, 1], p[, 2], lwd = lwd)
            }
            p <- rbind(to_xy(rc, hy(q[[1]])), to_xy(rc, hy(q[[2]])))
            graphics::lines(p[, 1], p[, 2], lwd = 0.8)
            p <- rbind(to_xy(rc, hy(q[[4]])), to_xy(rc, hy(q[[5]])))
            graphics::lines(p[, 1], p[, 2], lwd = 0.8)
            quad(rc - w * 0.22, rc + w * 0.22, hy(q[[2]]), hy(q[[4]]),
                border = "black", lwd = 0.9, col = NA)
            seg(rc - w * 0.22, rc + w * 0.22, hy(q[[3]]), lwd = 1.8)
            if (isTRUE(config$quartile_labels)) {
                labc <- c(pref("Diverging.RdBu:1"), pref("Diverging.RdBu:9"),
                    pref("Qualitative.Bold:2"))
                for (j in 2:4) {
                    xy <- to_xy(rc, hy(q[[j]]))
                    lab <- format(round(q[[j]], 1), nsmall = 1)
                    bw <- graphics::strwidth(lab, cex = 0.55) * 0.6
                    bh <- graphics::strheight(lab, cex = 0.55) * 0.75
                    graphics::rect(xy[1] - bw, xy[2] - bh, xy[1] + bw,
                        xy[2] + bh, col = "white", border = labc[[j - 1]],
                        lwd = 0.6)
                    graphics::text(xy[1], xy[2], lab, cex = 0.55,
                        col = labc[[j - 1]])
                }
            }
        } else {
            tab <- table(factor(t$v, levels = t$lv))
            pct <- as.numeric(tab) / sum(tab)
            h0 <- off + c(0, cumsum(pct))[seq_along(pct)] * H
            h1 <- off + cumsum(pct) * H
            for (j in seq_along(pct)) {
                if (pct[[j]] <= 0) next
                quad(r0, r1, h0[[j]], h1[[j]], col = t$pal[[j]],
                    border = "black", lwd = 0.5)
                if (pct[[j]] >= 0.05) {
                    xy <- to_xy((r0 + r1) / 2, (h0[[j]] + h1[[j]]) / 2)
                    graphics::text(xy[1], xy[2], sprintf("%.1f",
                        100 * pct[[j]]), cex = 0.5, col = "white",
                        srt = txt_srt)
                }
            }
        }
        tl <- to_xy((r0 + r1) / 2, off + H + 0.02)
        graphics::text(tl[1], tl[2], t$title, adj = c(0, 0.5), cex = 0.62,
            srt = txt_srt + 50)
    }

    tiles <- Filter(function(t) identical(t$type, "tile"), trk)
    if (length(tiles)) {
        lg <- lapply(tiles, function(t) {
            ComplexHeatmap::Legend(at = names(t$pal), nrow = 1,
                legend_gp = grid::gpar(fill = unname(t$pal)),
                title = if (is.null(t$legend)) t$column else t$legend,
                title_position = "leftcenter",
                title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
                labels_gp = grid::gpar(fontsize = 8),
                border = "black")
        })
        pk <- do.call(ComplexHeatmap::packLegend,
            c(lg, list(direction = "vertical", gap = grid::unit(2, "mm"))))
        ComplexHeatmap::draw(pk, x = grid::unit(config$legend_xy[[1]], "npc"),
            y = grid::unit(config$legend_xy[[2]], "npc"),
            just = c("left", "top"))
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
wd <- config$size$width
ht <- config$size$height
if (identical(fmt, "png")) {
    grDevices::png(out_path, width = wd, height = ht, units = "in",
        res = 300, bg = "transparent")
} else if (identical(fmt, "pdf")) {
    grDevices::pdf(out_path, width = wd, height = ht)
} else {
    grDevices::svg(out_path, width = wd, height = ht, bg = "transparent")
}
graphics::par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "transparent")
tryCatch(draw_figure(), finally = grDevices::dev.off())
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
