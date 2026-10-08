#!/usr/bin/env Rscript

# Template-ID: pie-doughnut-bar
#
# Purpose:
#   Draw a circos stacked doughnut of within-group composition on the
#   outer track, with radial bars of the same parts on the inner track.
#   One sector per group.
#
# Inputs:
#   One row per part within a group. Default example:
#     - group: sector (tissue, condition, cohort)
#     - item: part name within the group
#     - value: magnitude; turned into a within-group proportion
#
# Output:
#   A PDF, PNG, or SVG circular composition + bar figure.
#
# Dependencies:
#   circlize, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, inner_label, bar_height,
#   item_labels, sector_width, baseline). Equal vs count-weighted sector
#   width is CONFIG, not a second template. Item labels at the bar feet
#   are also CONFIG. A per-sample heat ring (SPI-style) is not this
#   template.
#
# Scientific assumptions:
#   value is a supplied magnitude, not computed in this script.
#   Doughnut slice angle is value / sum(value) within each group.
#   Bar height is the supplied value (default) or that within-group
#   proportion. Sorting by value inside a group is a display order,
#   not a statistical ranking.

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
        group = "group",
        item = "item",
        value = "value"
    ),
    inner_label = "value",
    bar_height = "value",
    item_labels = TRUE,
    sector_width = "count",
    baseline = 1,
    sort = "file",
    min_label = 0.05,
    size = list(
        width = 8,
        height = 8
    )
)

load_packages(c("circlize", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

g_col <- config$columns$group
i_col <- config$columns$item
v_col <- config$columns$value
df[[g_col]] <- as.character(df[[g_col]])
df[[i_col]] <- as.character(df[[i_col]])
df[[v_col]] <- as.numeric(df[[v_col]])
if (any(!is.finite(df[[v_col]]) | df[[v_col]] < 0)) {
    stop("value must be finite and non-negative.", call. = FALSE)
}

groups <- unique(df[[g_col]])
df[[g_col]] <- factor(df[[g_col]], levels = groups)

if (identical(config$sort, "value")) {
    df <- df[order(df[[g_col]], df[[v_col]]), , drop = FALSE]
} else {
    df <- df[order(df[[g_col]]), , drop = FALSE]
}

df$pct <- unlist(lapply(split(df[[v_col]], df[[g_col]]), function(v) {
    s <- sum(v)
    if (s <= 0) {
        stop("Each group needs a positive value sum.", call. = FALSE)
    }
    v / s
}), use.names = FALSE)
df$cumpct <- unlist(lapply(split(df$pct, df[[g_col]]), cumsum), use.names = FALSE)
df$n <- as.integer(table(df[[g_col]])[as.character(df[[g_col]])])
df$pos <- unlist(lapply(split(seq_len(nrow(df)), df[[g_col]]), function(idx) {
    n <- length(idx)
    (seq_len(n) - 0.5) / n
}), use.names = FALSE)

n_group <- length(groups)
# Lights from the supplied figure; darks are Brand.Elastic Search.
pair_light <- c(
    "#e3b5cb", "#fdf3c2", "#c5d69b",
    "#9cc1be", "#bed7eb", "#D5E8F1"
)
pair_dark <- palette_colors("Brand.Elastic Search")
n_pair <- min(length(pair_light), length(pair_dark))
color_maps <- vector("list", n_group)
names(color_maps) <- groups
for (i in seq_len(n_group)) {
    j <- ((i - 1L) %% n_pair) + 1L
    color_maps[[i]] <- c(pair_light[[j]], pair_dark[[j]])
}

sector_width <- if (identical(config$sector_width, "count")) {
    as.numeric(table(df[[g_col]])[groups])
} else {
    rep(1, n_group)
}

bar_y <- if (identical(config$bar_height, "percent")) {
    df$pct
} else {
    df[[v_col]]
}
max_y <- max(bar_y)
y_top <- max_y * 1.15
baseline <- config$baseline
if (!is.null(baseline) && (!is.finite(baseline) || baseline <= 0)) {
    baseline <- NULL
}
y_at <- if (!is.null(baseline) && baseline < y_top) {
    baseline
} else {
    at <- pretty(c(0, max_y), n = 4L)
    at[at > 0 & at < y_top]
}

inner_label <- config$inner_label
item_labels <- isTRUE(config$item_labels)
min_label <- as.numeric(config$min_label)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_circos <- function() {
    circlize::circos.clear()
    circlize::circos.par(
        start.degree = 90,
        gap.after = 1,
        track.margin = c(0, 0),
        points.overflow.warning = FALSE
    )
    circlize::circos.initialize(
        factors = groups,
        xlim = c(0, 1),
        sector.width = sector_width
    )

    circlize::circos.track(
        sectors = df[[g_col]],
        x = df$cumpct,
        ylim = c(0, 1),
        bg.border = NA,
        bg.col = NA,
        track.height = 0.22,
        panel.fun = function(x, y) {
            sector_name <- circlize::get.cell.meta.data("sector.index")
            pair <- color_maps[[sector_name]]
            ramp <- circlize::colorRamp2(c(0, 1), pair)
            left <- c(0, x[-length(x)])
            circlize::circos.rect(
                xleft = left,
                ybottom = 0,
                xright = x,
                ytop = 1,
                border = "white",
                col = ramp(x)
            )
            widths <- x - left
            labs <- rep(NA_character_, length(x))
            keep <- widths >= min_label
            if (identical(inner_label, "percent")) {
                labs[keep] <- sprintf("%.1f", 100 * widths[keep])
            } else if (identical(inner_label, "value")) {
                sub <- df[df[[g_col]] == sector_name, , drop = FALSE]
                labs[keep] <- sprintf("%s", round(sub[[v_col]][keep], 2))
            }
            if (identical(inner_label, "none")) {
                return(invisible())
            }
            circlize::circos.text(
                x = (left + x) / 2,
                y = 0.5,
                labels = labs,
                cex = 0.75,
                facing = "inside",
                niceFacing = TRUE
            )
        }
    )

    circlize::circos.track(
        sectors = df[[g_col]],
        x = df$pos,
        y = bar_y,
        ylim = c(0, y_top),
        track.height = 0.48,
        bg.border = NA,
        bg.col = NA,
        panel.fun = function(x, y) {
            sector_name <- circlize::get.cell.meta.data("sector.index")
            pair <- color_maps[[sector_name]]
            ramp <- circlize::colorRamp2(c(0, max_y), pair)
            n <- length(x)
            circlize::circos.barplot(
                y,
                x,
                bar_width = 0.72 / n,
                col = ramp(y)
            )
            if (!is.null(baseline) && baseline < y_top) {
                circlize::circos.lines(
                    x = c(0, 1),
                    y = c(baseline, baseline),
                    lwd = 0.8,
                    col = "black"
                )
            }
            if (length(y_at)) {
                circlize::circos.yaxis(
                    side = "right",
                    at = y_at,
                    labels.cex = 0.6,
                    labels.niceFacing = FALSE
                )
            }
            circlize::circos.text(
                x = 0.5,
                y = y_top * 0.92,
                labels = sector_name,
                facing = "bending.inside",
                niceFacing = TRUE,
                adj = c(0.5, 0.5),
                col = pair[[2]],
                cex = 1.1
            )
            if (item_labels) {
                sub <- df[df[[g_col]] == sector_name, , drop = FALSE]
                circlize::circos.text(
                    x = x,
                    y = -0.08 * y_top,
                    labels = sub[[i_col]],
                    adj = c(1, 0.5),
                    niceFacing = TRUE,
                    facing = "clockwise",
                    col = pair[[2]],
                    cex = 0.6
                )
            }
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
old_par <- graphics::par(mar = c(0.4, 0.4, 0.4, 0.4), bg = "transparent")
tryCatch(
    draw_circos(),
    finally = {
        graphics::par(old_par)
        grDevices::dev.off()
    }
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
