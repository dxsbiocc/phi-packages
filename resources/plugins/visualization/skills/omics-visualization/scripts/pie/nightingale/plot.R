#!/usr/bin/env Rscript

# Template-ID: pie-nightingale
#
# Purpose:
#   Draw a Nightingale rose: equal-angle polar bars, radius encodes
#   value, with rounded caps and an inner hole.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - roseName: category label
#     - value: numeric value
#
# Output:
#   A PDF, PNG, or SVG Nightingale rose chart.
#
# Dependencies:
#   ggplot2, readr, ggprism, gground, ggrepel, grid
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, inner hole, corner
#   shape = "petal" (example_petal.tsv) is the petal-pie variant:
#   equal-angle wedges from the centre, radius = value / max, each
#   ending in a semicircular cap, names along the radius, values at
#   the tip. Not a new id.
#   radius, value labels). Rounded caps vs a sharp rose, and the
#   inner hole, stay on this id. That is not pie-doughnut-round
#   (angle encodes share of a whole). Edit DATA PREPARATION to
#   change category order. Edit PLOT only when the glyph must change.
#
# Scientific assumptions:
#   Each row is one category with a comparable numeric value.
#   Equal angular width with radius encoding value is a display
#   choice. The inner hole and corner radius are display offsets,
#   not missing data. Values are already summarized.

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
        x = "roseName",
        y = "value"
    ),
    bar_width = 0.95,
    corner_pt = 5,
    inner_frac = 0.28,
    label_frac = 0.28,
    show_value_labels = TRUE,
    # "rounded" (default polar bars with rounded corners) or "petal":
    # equal-angle wedges meeting at the centre, radius by value, each
    # ending in a semicircular cap; names along the radius, values at
    # the tip, no legend.
    shape = "rounded",
    petal = list(
        palette = "Qualitative.Vivid / 10",
        min_frac = 0.2,          # floor for wedge length (radius = value / max)
        border = "white",
        name_size = 4.2,
        value_size = 4.4,
        start = 90               # degrees; first row starts here, clockwise
    ),
    labels = list(
        title = "",
        x = "",
        y = "",
        fill = NULL
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "gground", "ggrepel"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])
ok <- is.finite(df[[y_col]]) & df[[y_col]] >= 0
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with invalid values.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

y_max <- max(df[[y_col]], na.rm = TRUE)
if (!is.finite(y_max) || y_max <= 0) {
    stop("Need a positive numeric value for the rose radius.", call. = FALSE)
}
inner_frac <- as.numeric(config$inner_frac)
if (!is.finite(inner_frac) || inner_frac < 0) {
    inner_frac <- 0.28
}
label_frac <- as.numeric(config$label_frac)
if (!is.finite(label_frac) || label_frac < 0) {
    label_frac <- 0.28
}
inner <- y_max * inner_frac
y_top <- y_max + y_max * label_frac
bar_width <- as.numeric(config$bar_width)
if (!is.finite(bar_width) || bar_width <= 0 || bar_width > 1) {
    bar_width <- 0.95
}
corner_pt <- as.numeric(config$corner_pt)
if (!is.finite(corner_pt) || corner_pt < 0) {
    corner_pt <- 5
}

# Reference rose fills from the source figure (圆角扇形图).
rose_seq <- c(
    "#5070dd", "#b6d634", "#ff994d", "#0ca8df",
    "#ffd10a", "#fb628b", "#785db0", "#3fbe95"
)
n_fill <- nlevels(df[[x_col]])
if (n_fill <= length(rose_seq)) {
    fill_vals <- rose_seq[seq_len(n_fill)]
} else {
    fill_vals <- palette_colors("Brand.Algolia", n = n_fill)
    fill_vals[seq_along(rose_seq)] <- rose_seq
}
ink <- palette_colors("Brand.Emma")[[1]]
seg <- palette_colors("Qualitative.Safe")[[12]]

if (identical(config$shape, "petal")) {
    pc <- config$petal
    n <- nrow(df)
    v <- df[[y_col]]
    r_side <- pmax(v / max(v), pc$min_frac)
    slot <- 2 * pi / n
    a0 <- pc$start / 180 * pi - (seq_len(n) - 1) * slot
    a1 <- a0 - slot
    polys <- do.call(rbind, lapply(seq_len(n), function(i) {
        r <- r_side[[i]]
        p0 <- r * c(cos(a0[[i]]), sin(a0[[i]]))
        p1 <- r * c(cos(a1[[i]]), sin(a1[[i]]))
        m <- (p0 + p1) / 2
        rc <- sqrt(sum((p0 - m)^2))
        mid <- (a0[[i]] + a1[[i]]) / 2
        ang0 <- atan2(p0[[2]] - m[[2]], p0[[1]] - m[[1]])
        cap <- seq(ang0, ang0 - pi, length.out = 60)
        cx <- m[[1]] + rc * cos(cap)
        cy <- m[[2]] + rc * sin(cap)
        # keep the cap bulging outward
        if (sum((c(mean(cx), mean(cy)) - m) * c(cos(mid), sin(mid))) < 0) {
            cap <- seq(ang0, ang0 + pi, length.out = 60)
            cx <- m[[1]] + rc * cos(cap)
            cy <- m[[2]] + rc * sin(cap)
        }
        data.frame(id = i, x = c(0, cx, 0), y = c(0, cy, 0))
    }))
    mid <- (a0 + a1) / 2
    deg <- mid * 180 / pi
    left <- cos(mid) < -1e-9
    reach <- r_side * cos(slot / 2) + r_side * sin(slot / 2)
    lab <- data.frame(
        id = seq_len(n),
        name = as.character(df[[x_col]]), value = v,
        nx = 0.55 * r_side * cos(mid), ny = 0.55 * r_side * sin(mid),
        vx = (reach - 0.12) * cos(mid), vy = (reach - 0.12) * sin(mid),
        nrot = ifelse(left, deg + 180, deg),
        vrot = deg - 90
    )
    cols <- palette_colors(pc$palette, n = n)
    p <- ggplot() +
        geom_polygon(data = polys, aes(x = x, y = y, group = id,
            fill = factor(id)), colour = pc$border, linewidth = 0.6) +
        geom_text(data = lab, aes(x = nx, y = ny, label = name, angle = nrot),
            colour = "white", size = pc$name_size) +
        scale_fill_manual(values = cols, guide = "none") +
        coord_fixed(clip = "off") +
        labs(title = config$labels$title) +
        theme_void() +
        theme(plot.background = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    if (isTRUE(config$show_value_labels)) {
        p <- p + geom_text(data = lab, aes(x = vx, y = vy, label = value,
            angle = vrot), colour = "white", size = pc$value_size)
    }
    save_ggplot(p, io$output, width = 7, height = 7)
    quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[x_col]]
)) +
    geom_round_col(
        width = bar_width,
        radius = grid::unit(corner_pt, "pt"),
        colour = NA
    ) +
    scale_fill_manual(name = config$labels$fill, values = fill_vals) +
    coord_polar(clip = "off") +
    scale_y_continuous(
        limits = c(-inner, y_top),
        expand = c(0, 0)
    ) +
    labs(
        title = config$labels$title,
        x = NULL,
        y = NULL
    ) +
    guides(fill = guide_legend(nrow = 1, title = NULL)) +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.position = "bottom",
        legend.key = element_blank()
    )

if (isTRUE(config$show_value_labels)) {
    p <- p + geom_text_repel(
        aes(label = .data[[y_col]]),
        size = 3.6,
        nudge_y = y_max * 0.22,
        segment.curvature = -0.1,
        segment.ncp = 3,
        segment.angle = 20,
        segment.size = 0.4,
        segment.colour = seg,
        colour = ink,
        min.segment.length = 0,
        show.legend = FALSE
    )
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 7.2)
