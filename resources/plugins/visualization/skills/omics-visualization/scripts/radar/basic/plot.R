#!/usr/bin/env Rscript

# Template-ID: radar-basic
#
# Purpose:
#   Draw a radar chart comparing a few series across named axes.
#
# Inputs:
#   A wide table with one row per axis. Default example:
#     - name: axis label ("\n" breaks a line)
#     - remaining columns: one numeric series each (Budget, Spending)
#     - optional axis_max / axis_group columns (see CONFIG); they are
#       not drawn as series
#
# Output:
#   A PDF, PNG, or SVG radar chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#   example_units.tsv needs axis_scale = "per_axis" in CONFIG.
#   example_models.tsv needs r_min, shapes and a per-axis band in CONFIG.
#
# Agent adaptation:
#   For a new table, edit only CONFIG. Switches:
#     axis_scale    "shared" (one radial scale) or "per_axis" (each axis
#                   scaled to its own max; use when units differ)
#     r_min         radial baseline for shared scale; > 0 truncates the
#                   centre and forces radial tick labels
#     series_shapes one point shape per series
#     smooth        periodic-spline outline (display only); show_points
#     log_axes      axis names drawn on log10 before per-axis scaling
#     value_labels  print raw values next to vertices
#     grid_shape    "polygon" or "circle"; grid_linetype e.g. "dotted"
#     axis_band     "auto" (axis_group column), "per_axis" (one coloured
#                   segment per axis holding its name), or "none"
#     direction / start_degree   axis order and rotation
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one axis; remaining columns are series.
#   With axis_scale = "shared", radius is comparable across axes. With
#   r_min > 0 the centre is truncated: distance and area are no longer
#   proportional to value, so radial ticks are always printed.
#   With axis_scale = "per_axis", radius is comparable only along one
#   axis; read magnitudes from the printed raw values. log_axes is a
#   display transform. Closing the polygon repeats the first axis.
#   smooth = TRUE draws a periodic spline through the vertices, clamped at
#   0; axes are discrete, so the curve between two axes is not data and
#   may overshoot. Points mark the real values.
#   Values are already summarized, not raw replicates.

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
        axis = "name",
        axis_max = "axis_max",
        axis_group = "axis_group"
    ),
    labels = list(
        title = "Radar Chart",
        axis_group = ""
    ),
    series_colours = c("#037ef3", "#f85a40"),
    # Point shape per series (e.g. c(17, 8, 15, 18, 16)); NULL = filled circles.
    series_shapes = NULL,
    # Periodic-spline outline through the vertices (display only; the
    # curve between axes is not data). show_points keeps the real values.
    smooth = FALSE,
    show_points = TRUE,
    axis_scale = "shared",
    # Radial baseline for axis_scale = "shared". > 0 truncates the
    # centre, so radial tick values are always printed.
    r_min = 0,
    tick_labels = "auto",       # "auto" (only when r_min > 0), TRUE, FALSE
    tick_angle_deg = NULL,      # NULL = halfway between the first two axes
    log_axes = character(0),
    value_labels = FALSE,
    value_digits = 3,
    grid_shape = "polygon",
    grid_linetype = "solid",
    n_rings = 5L,
    ring_fills = c("#E8F8F2", "#FDF6E6"),
    fill_alpha = 0.2,
    direction = "counterclockwise",
    start_degree = 0,           # rotate the first axis away from 12 o'clock
    axis_label_wrap = 18L,
    # "auto" = band by axis_group column if present; "per_axis" = one
    # coloured segment per axis with the axis name inside; "none".
    axis_band = "auto",
    axis_group_fills = c("#cde8ef", "#e7e5e8", "#f6e3c9", "#dcefd8"),
    axis_band_fills = c("#95BBE1", "#C6BAE2", "#EE8D81", "#FC9D4F",
        "#F9E0CC", "#529BB3", "#F9DA7B"),
    axis_label_size = 3.2,
    band_label_size = 3.6,
    value_label_size = 3.0,
    tick_label_size = 3.0
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns["axis"])
axis_col <- config$columns$axis
max_col <- config$columns$axis_max
grp_col <- config$columns$axis_group
has_max <- !is.null(max_col) && max_col %in% names(df)
has_grp <- !is.null(grp_col) && grp_col %in% names(df)
series_cols <- setdiff(names(df), c(axis_col, if (has_max) max_col, if (has_grp) grp_col))
if (!length(series_cols)) {
    stop("No series columns found.", call. = FALSE)
}

axis_names <- gsub("\\\\n", "\n", as.character(df[[axis_col]]))
n_axes <- length(axis_names)
if (n_axes < 3L) {
    stop("A radar needs at least three axes.", call. = FALSE)
}
vals <- as.matrix(df[series_cols])
storage.mode(vals) <- "numeric"

axis_scale <- match.arg(config$axis_scale, c("shared", "per_axis"))
grid_shape <- match.arg(config$grid_shape, c("polygon", "circle"))
is_log <- axis_names %in% config$log_axes | as.character(df[[axis_col]]) %in% config$log_axes

axis_band <- match.arg(config$axis_band, c("auto", "per_axis", "none"))
if (identical(axis_band, "auto")) {
    axis_band <- if (has_grp) "group" else "none"
}
r_min <- if (identical(axis_scale, "shared")) as.numeric(config$r_min) else 0
breaks <- NULL

# Radius in [0, 1] for every vertex.
if (identical(axis_scale, "shared")) {
    v_hi <- max(vals, na.rm = TRUE)
    if (r_min > 0 && any(vals < r_min, na.rm = TRUE)) {
        warning("Some values are below r_min and are drawn at the centre.",
            call. = FALSE)
    }
    breaks <- pretty(c(r_min, v_hi), n = config$n_rings)
    breaks <- breaks[breaks > r_min]
    y_max <- max(breaks)
    if (!is.finite(y_max) || y_max <= r_min) {
        y_max <- r_min + 1
        breaks <- y_max
    }
    rad <- (vals - r_min) / (y_max - r_min)
} else {
    tv <- vals
    if (any(is_log)) {
        if (any(vals[is_log, ] <= 0, na.rm = TRUE)) {
            stop("log_axes need positive values.", call. = FALSE)
        }
        tv[is_log, ] <- log10(vals[is_log, ])
    }
    amax <- apply(tv, 1, max, na.rm = TRUE)
    if (has_max) {
        sup <- as.numeric(df[[max_col]])
        sup[is_log] <- log10(sup[is_log])
        amax <- ifelse(is.finite(sup) & sup > 0, sup, amax)
    }
    amax[!is.finite(amax) | amax <= 0] <- 1
    rad <- tv / amax
}
rad[rad < 0] <- 0

dir_sign <- if (identical(config$direction, "clockwise")) -1 else 1
theta <- pi / 2 - config$start_degree * pi / 180 +
    dir_sign * 2 * pi * (seq_len(n_axes) - 1) / n_axes
long <- data.frame(
    idx = rep(seq_len(n_axes), times = length(series_cols)),
    series = factor(rep(series_cols, each = n_axes), levels = series_cols),
    value = as.numeric(vals),
    r = as.numeric(rad),
    stringsAsFactors = FALSE
)
long$px <- long$r * cos(theta[long$idx])
long$py <- long$r * sin(theta[long$idx])
long_closed <- rbind(long, long[long$idx == 1L, , drop = FALSE])

# Smooth outline: periodic spline of radius over axis position, clamped at
# 0 so it never dips below the centre. Display only; points stay on data.
if (isTRUE(config$smooth)) {
    n_fine <- n_axes * 24L
    long_closed <- do.call(rbind, lapply(series_cols, function(s) {
        r <- long$r[long$series == s]
        sp <- stats::spline(seq_len(n_axes + 1L), c(r, r[[1]]),
            n = n_fine + 1L, method = "periodic")
        a <- pi / 2 - config$start_degree * pi / 180 +
            dir_sign * 2 * pi * (sp$x - 1) / n_axes
        rs <- pmax(sp$y, 0)
        data.frame(series = factor(s, levels = series_cols),
            px = rs * cos(a), py = rs * sin(a))
    }))
}

# Grid rings.
ring_path <- function(r) {
    if (identical(grid_shape, "circle")) {
        a <- seq(0, 2 * pi, length.out = 181)
    } else {
        a <- c(theta, theta[[1]])
    }
    data.frame(px = r * cos(a), py = r * sin(a))
}
ring_r <- if (is.null(breaks)) {
    seq_len(config$n_rings) / config$n_rings
} else {
    (breaks - r_min) / (y_max - r_min)
}
n_ring <- length(ring_r)
rings <- do.call(rbind, lapply(seq_len(n_ring), function(i) {
    d <- ring_path(ring_r[[i]])
    d$ring <- i
    d
}))
ring_geoms <- lapply(rev(seq_len(n_ring)), function(i) {
    ggplot2::geom_polygon(
        data = rings[rings$ring == i, , drop = FALSE],
        ggplot2::aes(x = px, y = py),
        fill = config$ring_fills[[((i - 1L) %% length(config$ring_fills)) + 1L]],
        colour = NA,
        inherit.aes = FALSE
    )
})
spokes <- data.frame(px = cos(theta), py = sin(theta))

# Radial tick values (forced on when the centre is truncated).
show_ticks <- if (identical(config$tick_labels, "auto")) {
    r_min > 0
} else {
    isTRUE(config$tick_labels)
}
tick_df <- NULL
if (show_ticks && !is.null(breaks)) {
    step0 <- 2 * pi / n_axes
    ta <- if (is.null(config$tick_angle_deg)) {
        theta[[1]] + dir_sign * step0 / 2
    } else {
        config$tick_angle_deg * pi / 180
    }
    tick_df <- data.frame(
        px = ring_r * cos(ta),
        py = ring_r * sin(ta),
        label = format(breaks, trim = TRUE, drop0trailing = FALSE),
        stringsAsFactors = FALSE
    )
}

# Axis labels: wrapped, anchored by side so nothing is clipped.
wrap_lab <- function(s) {
    vapply(strsplit(s, "\n", fixed = TRUE), function(parts) {
        paste(unlist(lapply(parts, strwrap, width = config$axis_label_wrap)),
            collapse = "\n")
    }, character(1))
}
band_in <- 1.04
band_out <- if (identical(axis_band, "per_axis")) 1.18 else 1.12
lab_r <- (if (identical(axis_band, "none")) 1.10 else band_out + 0.08) +
    (if (isTRUE(config$value_labels)) 0.10 else 0)
axis_lab <- data.frame(
    px = lab_r * cos(theta),
    py = lab_r * sin(theta),
    label = wrap_lab(axis_names),
    stringsAsFactors = FALSE
)
axis_lab$hjust <- ifelse(abs(axis_lab$px) < 0.08, 0.5, ifelse(axis_lab$px < 0, 1, 0))
axis_lab$vjust <- ifelse(abs(axis_lab$py) < 0.08, 0.5, ifelse(axis_lab$py < 0, 1, 0))

# Raw value labels: outermost series per axis outside, the rest inside.
fmt_value <- function(x) {
    sup <- c("0" = "\u2070", "1" = "\u00b9", "2" = "\u00b2", "3" = "\u00b3",
        "4" = "\u2074", "5" = "\u2075", "6" = "\u2076", "7" = "\u2077",
        "8" = "\u2078", "9" = "\u2079", "-" = "\u207b")
    vapply(x, function(v) {
        if (!is.finite(v)) return("")
        if (abs(v) >= 1e5 || (abs(v) > 0 && abs(v) < 1e-3)) {
            e <- floor(log10(abs(v)))
            m <- signif(v / 10^e, config$value_digits)
            ex <- paste(sup[strsplit(as.character(e), "")[[1]]], collapse = "")
            if (isTRUE(all.equal(m, 1))) paste0("10", ex) else paste0(m, "\u00d710", ex)
        } else {
            formatC(signif(v, config$value_digits), format = "fg", big.mark = ",",
                digits = config$value_digits)
        }
    }, character(1))
}
val_lab <- NULL
if (isTRUE(config$value_labels)) {
    val_lab <- long
    top <- stats::ave(val_lab$r, val_lab$idx, FUN = max)
    outside <- val_lab$r >= top - 1e-9
    rr <- ifelse(outside, val_lab$r + 0.08, pmax(val_lab$r - 0.10, 0.05))
    val_lab$px <- rr * cos(theta[val_lab$idx])
    val_lab$py <- rr * sin(theta[val_lab$idx])
    val_lab$label <- trimws(fmt_value(val_lab$value))
}

# Outer band: by axis_group column, or one segment per axis.
tangent_angle <- function(th) {
    a <- th * 180 / pi - 90
    a <- ((a + 180) %% 360) - 180
    ifelse(a > 90 | a < -90, a + 180, a)
}
band <- NULL
if (!identical(axis_band, "none")) {
    step <- 2 * pi / n_axes
    if (identical(axis_band, "per_axis")) {
        g <- paste0("axis", seq_len(n_axes))
        half <- step / 2 * 0.97
    } else {
        g <- as.character(df[[grp_col]])
        half <- step / 2
    }
    g_levels <- unique(g)
    band <- do.call(rbind, lapply(seq_len(n_axes), function(i) {
        a <- seq(theta[[i]] - half, theta[[i]] + half, length.out = 30)
        data.frame(
            px = c(band_in * cos(a), band_out * cos(rev(a))),
            py = c(band_in * sin(a), band_out * sin(rev(a))),
            seg = i,
            axis_group = factor(g[[i]], levels = g_levels),
            stringsAsFactors = FALSE
        )
    }))
    band_pal <- if (identical(axis_band, "per_axis")) {
        config$axis_band_fills
    } else {
        config$axis_group_fills
    }
    band_fills <- expand_palette(band_pal, length(g_levels))
    names(band_fills) <- g_levels
    if (identical(axis_band, "per_axis")) {
        band_lab <- data.frame(
            px = (band_in + band_out) / 2 * cos(theta),
            py = (band_in + band_out) / 2 * sin(theta),
            label = gsub("\n", " ", axis_names),
            angle = tangent_angle(theta),
            stringsAsFactors = FALSE
        )
    } else {
        # One text label per contiguous run of the same group.
        runs <- rle(g)
        ends <- cumsum(runs$lengths)
        starts <- ends - runs$lengths + 1L
        mid_theta <- vapply(seq_along(starts), function(k) {
            mean(theta[starts[[k]]:ends[[k]]])
        }, numeric(1))
        band_lab <- data.frame(
            px = (band_in + band_out) / 2 * cos(mid_theta),
            py = (band_in + band_out) / 2 * sin(mid_theta),
            label = runs$values,
            angle = tangent_angle(mid_theta),
            stringsAsFactors = FALSE
        )
    }
}

series_cols_hex <- expand_palette(config$series_colours, length(series_cols))
names(series_cols_hex) <- series_cols
shape_vals <- if (is.null(config$series_shapes)) {
    rep(16, length(series_cols))
} else {
    rep_len(config$series_shapes, length(series_cols))
}
names(shape_vals) <- series_cols

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot2::ggplot()
if (!is.null(band)) {
    for (lv in names(band_fills)) {
        p <- p + ggplot2::geom_polygon(
            data = band[band$axis_group == lv, , drop = FALSE],
            ggplot2::aes(x = px, y = py, group = seg),
            fill = band_fills[[lv]],
            colour = NA
        )
    }
}
p <- p +
    ring_geoms +
    ggplot2::geom_path(
        data = rings,
        ggplot2::aes(x = px, y = py, group = ring),
        colour = if (identical(config$grid_linetype, "solid")) "#dbdee4" else "grey45",
        linewidth = 0.35,
        linetype = config$grid_linetype
    ) +
    ggplot2::geom_segment(
        data = spokes,
        ggplot2::aes(x = 0, y = 0, xend = px, yend = py),
        colour = "#cfd2d7",
        linewidth = 0.4
    ) +
    (if (config$fill_alpha > 0) {
        ggplot2::geom_polygon(
            data = long_closed,
            ggplot2::aes(x = px, y = py, colour = series, fill = series, group = series),
            alpha = config$fill_alpha,
            linewidth = 0.8
        )
    } else {
        ggplot2::geom_path(
            data = long_closed,
            ggplot2::aes(x = px, y = py, colour = series, group = series),
            linewidth = 0.8,
            linejoin = "round"
        )
    }) +
    ggplot2::geom_point(
        data = long,
        ggplot2::aes(x = px, y = py, colour = series, shape = series),
        size = if (is.null(config$series_shapes)) 2 else 2.8,
        alpha = if (isTRUE(config$show_points)) 1 else 0
    ) +
    ggplot2::scale_colour_manual(values = series_cols_hex, name = NULL) +
    ggplot2::scale_fill_manual(
        values = series_cols_hex, name = NULL,
        guide = if (config$fill_alpha > 0) "legend" else "none"
    ) +
    ggplot2::scale_shape_manual(values = shape_vals, name = NULL)

if (!identical(axis_band, "per_axis")) {
    p <- p + ggplot2::geom_text(
        data = axis_lab,
        ggplot2::aes(x = px, y = py, label = label, hjust = hjust, vjust = vjust),
        size = config$axis_label_size,
        lineheight = 0.9
    )
}

if (!is.null(band)) {
    per_axis_band <- identical(axis_band, "per_axis")
    p <- p + ggplot2::geom_text(
        data = band_lab,
        ggplot2::aes(x = px, y = py, label = label, angle = angle),
        size = if (per_axis_band) config$band_label_size else 2.6,
        colour = if (per_axis_band) "grey10" else "grey30",
        fontface = if (per_axis_band) "plain" else "italic"
    )
}

if (!is.null(tick_df)) {
    p <- p + ggplot2::geom_label(
        data = tick_df,
        ggplot2::aes(x = px, y = py, label = label),
        size = config$tick_label_size,
        fontface = "bold",
        fill = "white",
        label.padding = ggplot2::unit(0.08, "lines"),
        border.colour = NA
    )
}

if (!is.null(val_lab)) {
    p <- p + ggplot2::geom_text(
        data = val_lab,
        ggplot2::aes(x = px, y = py, label = label, colour = series),
        size = config$value_label_size,
        fontface = "bold",
        show.legend = FALSE
    )
}

lim <- if (identical(axis_band, "per_axis")) band_out + 0.04 else lab_r + 0.05
p <- p +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim),
        expand = FALSE, clip = "off") +
    ggplot2::labs(title = config$labels$title, x = NULL, y = NULL) +
    theme_prism() +
    ggplot2::theme(
        plot.margin = ggplot2::margin(20, 20, 20, 90),
        # Axis labels overhang the panel (clip = "off"); keep the legend clear.
        legend.box.spacing = ggplot2::unit(30, "mm"),
        plot.background = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, margin = ggplot2::margin(b = 24)),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 9, height = 8)
