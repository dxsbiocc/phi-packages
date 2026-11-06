#!/usr/bin/env Rscript

# Template-ID: radar-group
#
# Purpose:
#   Draw one radar chart per series in a faceted grid (not overlaid).
#
# Inputs:
#   A wide table with one row per axis. Default example:
#     - name: axis label ("\n" breaks a line)
#     - remaining columns: one numeric series per panel
#     - optional axis_max / axis_group columns (see CONFIG); they are
#       not drawn as series
#
# Output:
#   A PDF, PNG, or SVG grouped radar chart.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggh4x
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG. Switches:
#     axis_scale   "shared" (one radial scale) or "per_axis" (each axis
#                  scaled to its own max; use when units differ)
#     log_axes     axis names drawn on log10 before per-axis scaling
#     value_labels print raw values next to vertices
#     grid_shape   "polygon" or "circle"
#     strip_style  "plain", "box" (filled strip), "text" (series-coloured title)
#     axis_group   column that groups axes; draws an outer colour band
#     fill_mode    "centre" (fill to the line) or "baseline" (fill the area
#                  between a baseline and the line, cut at crossings);
#                  baseline is a scalar in CONFIG or a per-axis column
#     point_style  "filled" or "hollow"; outer_frame darkens the outer ring
#     ncol         panels per row
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each series is drawn in its own panel. All panels share the same
#   axes and the same radial scaling, so panels stay comparable.
#   With axis_scale = "per_axis", radius is comparable only along one
#   axis; read magnitudes from the printed raw values. log_axes is a
#   display transform. Closing the polygon repeats the first axis.
#   With fill_mode = "baseline" the shaded area shows deviation from a
#   supplied baseline (e.g. chance level 1 / n axes); it is not a test.
#   Scores are already computed, not raw replicates.

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
        axis_group = "axis_group",
        baseline = "baseline"       # optional per-axis baseline column
    ),
    # "centre": fill from the centre to the line. "baseline": fill the area
    # between the baseline and the line (above and below).
    fill_mode = "auto",              # "auto" = baseline when one is supplied, else centre
    baseline = NULL,                # scalar baseline in data units; column wins
    baseline_line = TRUE,
    strip_style = "text",           # "plain", "box" (filled strip), "text" (series-coloured title)
    outer_frame = TRUE,             # dark outline on the outermost ring
    point_style = "hollow",         # "filled" or "hollow" (white centre)
    labels = list(
        title = "Importance of facial features",
        axis_group = ""
    ),
    series_colours = c("#1d6996", "#6699cc", "#888888", "#cc503e"),
    axis_scale = "shared",
    log_axes = character(0),
    value_labels = FALSE,
    value_digits = 3,
    grid_shape = "circle",
    strip_colour = FALSE,
    ncol = 3L,
    n_rings = 5L,
    ring_fills = c("white"),
    fill_alpha = 0.5,
    axis_label_wrap = 14L,
    axis_group_fills = c("#cde8ef", "#e7e5e8", "#f6e3c9", "#dcefd8"),
    axis_label_size = 2.9,
    value_label_size = 2.6
)

load_packages(c("ggplot2", "readr", "ggprism", "ggh4x"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns["axis"])
axis_col <- config$columns$axis
max_col <- config$columns$axis_max
grp_col <- config$columns$axis_group
base_col <- config$columns$baseline
has_max <- !is.null(max_col) && max_col %in% names(df)
has_grp <- !is.null(grp_col) && grp_col %in% names(df)
has_base_col <- !is.null(base_col) && base_col %in% names(df)
series_cols <- setdiff(names(df), c(axis_col, if (has_max) max_col, if (has_grp) grp_col,
    if (has_base_col) base_col))
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
fill_mode <- match.arg(config$fill_mode, c("auto", "centre", "baseline"))
is_log <- axis_names %in% config$log_axes | as.character(df[[axis_col]]) %in% config$log_axes

base_val <- if (has_base_col) {
    as.numeric(df[[base_col]])
} else if (!is.null(config$baseline)) {
    rep(as.numeric(config$baseline), n_axes)
} else {
    NULL
}
if (identical(fill_mode, "auto")) {
    fill_mode <- if (is.null(base_val)) "centre" else "baseline"
}
if (identical(fill_mode, "baseline") && is.null(base_val)) {
    stop("fill_mode = \"baseline\" needs config$baseline or a baseline column.", call. = FALSE)
}

# Radius in [0, 1] for every vertex.
if (identical(axis_scale, "shared")) {
    y_max <- max(pretty(c(0, max(c(vals, base_val), na.rm = TRUE))))
    if (!is.finite(y_max) || y_max <= 0) {
        y_max <- 1
    }
    rad <- vals / y_max
    base_r <- if (is.null(base_val)) NULL else base_val / y_max
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
    base_r <- if (is.null(base_val)) {
        NULL
    } else {
        bv <- base_val
        bv[is_log] <- log10(bv[is_log])
        bv / amax
    }
}
rad[rad < 0] <- 0
if (!is.null(base_r)) base_r <- pmax(base_r, 0)

theta <- pi / 2 + 2 * pi * (seq_len(n_axes) - 1) / n_axes
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

# Baseline curve: a circle when constant, else straight segments between
# per-axis baseline points. base_radius() gives its radius at angle phi
# inside sector i -> i + 1.
base_const <- !is.null(base_r) && diff(range(base_r)) < 1e-12
base_radius <- function(phi, i) {
    if (base_const) return(rep(base_r[[1]], length(phi)))
    j <- if (i == n_axes) 1L else i + 1L
    a <- c(base_r[[i]] * cos(theta[[i]]), base_r[[i]] * sin(theta[[i]]))
    b <- c(base_r[[j]] * cos(theta[[j]]), base_r[[j]] * sin(theta[[j]]))
    # Ray (cos phi, sin phi) * r meets segment a + s (b - a).
    d <- b - a
    den <- cos(phi) * d[[2]] - sin(phi) * d[[1]]
    (a[[1]] * d[[2]] - a[[2]] * d[[1]]) / den
}
base_path <- NULL
if (!is.null(base_r)) {
    base_path <- do.call(rbind, lapply(seq_len(n_axes), function(i) {
        j <- if (i == n_axes) 1L else i + 1L
        t2 <- theta[[j]] + if (j == 1L) 2 * pi else 0
        phi <- seq(theta[[i]], t2, length.out = 30)
        r <- base_radius(phi, i)
        data.frame(px = r * cos(phi), py = r * sin(phi))
    }))
}

# Band between the baseline and the line, cut at every crossing so above
# and below pieces never bridge the intersection.
band_pieces <- NULL
if (identical(fill_mode, "baseline")) {
    tt <- seq(0, 1, length.out = 41)
    band_pieces <- do.call(rbind, lapply(series_cols, function(s) {
        ls <- long[long$series == s, , drop = FALSE]
        do.call(rbind, lapply(seq_len(n_axes), function(i) {
            j <- if (i == n_axes) 1L else i + 1L
            x <- ls$px[[i]] + (ls$px[[j]] - ls$px[[i]]) * tt
            y <- ls$py[[i]] + (ls$py[[j]] - ls$py[[i]]) * tt
            phi <- atan2(y, x)
            # Keep phi continuous across the sector.
            ref <- theta[[i]]
            phi <- ref + ((phi - ref + pi) %% (2 * pi)) - pi
            rb <- base_radius(phi, i)
            sgn <- sqrt(x^2 + y^2) - rb
            # Split at sign changes (linear root between samples).
            cut <- which(diff(sign(sgn)) != 0)
            roots <- vapply(cut, function(c0) {
                f <- sgn[[c0]] / (sgn[[c0]] - sgn[[c0 + 1]])
                tt[[c0]] + f * (tt[[c0 + 1]] - tt[[c0]])
            }, numeric(1))
            edges <- c(0, roots, 1)
            pieces <- list()
            for (k in seq_len(length(edges) - 1L)) {
                start_t <- edges[[k]]
                end_t <- edges[[k + 1L]]
                seg_t <- unique(c(start_t, tt[tt > start_t & tt < end_t], end_t))
                sx <- ls$px[[i]] + (ls$px[[j]] - ls$px[[i]]) * seg_t
                sy <- ls$py[[i]] + (ls$py[[j]] - ls$py[[i]]) * seg_t
                sphi <- atan2(sy, sx)
                sphi <- ref + ((sphi - ref + pi) %% (2 * pi)) - pi
                srb <- base_radius(sphi, i)
                pieces[[k]] <- data.frame(
                    px = c(sx, rev(srb * cos(sphi))),
                    py = c(sy, rev(srb * sin(sphi))),
                    id = paste(s, i, k, sep = "_"),
                    series = s,
                    stringsAsFactors = FALSE
                )
            }
            do.call(rbind, pieces)
        }))
    }))
    band_pieces$series <- factor(band_pieces$series, levels = series_cols)
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
rings <- do.call(rbind, lapply(seq_len(config$n_rings), function(i) {
    d <- ring_path(i / config$n_rings)
    d$ring <- i
    d
}))
ring_geoms <- lapply(rev(seq_len(config$n_rings)), function(i) {
    ggplot2::geom_polygon(
        data = rings[rings$ring == i, , drop = FALSE],
        ggplot2::aes(x = px, y = py),
        fill = config$ring_fills[[((i - 1L) %% length(config$ring_fills)) + 1L]],
        colour = NA,
        inherit.aes = FALSE
    )
})
spokes <- data.frame(px = cos(theta), py = sin(theta))

# Axis labels: wrapped, anchored by side so nothing is clipped.
wrap_lab <- function(s) {
    vapply(strsplit(s, "\n", fixed = TRUE), function(parts) {
        paste(unlist(lapply(parts, strwrap, width = config$axis_label_wrap)),
            collapse = "\n")
    }, character(1))
}
lab_r <- (if (has_grp) 1.20 else 1.10) + (if (isTRUE(config$value_labels)) 0.10 else 0)
axis_lab <- data.frame(
    px = lab_r * cos(theta),
    py = lab_r * sin(theta),
    label = wrap_lab(axis_names),
    stringsAsFactors = FALSE
)
axis_lab$hjust <- ifelse(abs(axis_lab$px) < 0.08, 0.5, ifelse(axis_lab$px < 0, 1, 0))
axis_lab$vjust <- ifelse(abs(axis_lab$py) < 0.08, 0.5, ifelse(axis_lab$py < 0, 1, 0))

# Raw value labels, just outside each vertex (one series per panel).
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
    rr <- val_lab$r + 0.10
    val_lab$px <- rr * cos(theta[val_lab$idx])
    val_lab$py <- rr * sin(theta[val_lab$idx])
    val_lab$label <- trimws(fmt_value(val_lab$value))
}

# Outer axis-group band.
band <- NULL
if (has_grp) {
    g <- as.character(df[[grp_col]])
    g_levels <- unique(g)
    step <- 2 * pi / n_axes
    band <- do.call(rbind, lapply(seq_len(n_axes), function(i) {
        a <- seq(theta[[i]] - step / 2, theta[[i]] + step / 2, length.out = 30)
        data.frame(
            px = c(1.04 * cos(a), 1.12 * cos(rev(a))),
            py = c(1.04 * sin(a), 1.12 * sin(rev(a))),
            seg = i,
            axis_group = factor(g[[i]], levels = g_levels),
            stringsAsFactors = FALSE
        )
    }))
    band_fills <- expand_palette(config$axis_group_fills, length(g_levels))
    names(band_fills) <- g_levels
    # One text label per contiguous run of the same group, along the band.
    runs <- rle(g)
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L
    mid_theta <- vapply(seq_along(starts), function(k) {
        mean(theta[starts[[k]]:ends[[k]]])
    }, numeric(1))
    band_lab <- data.frame(
        px = 1.08 * cos(mid_theta),
        py = 1.08 * sin(mid_theta),
        label = runs$values,
        angle = {
            a <- mid_theta * 180 / pi - 90
            a <- ((a + 180) %% 360) - 180
            ifelse(a > 90 | a < -90, a + 180, a)
        },
        stringsAsFactors = FALSE
    )
}

series_cols_hex <- expand_palette(config$series_colours, length(series_cols))
names(series_cols_hex) <- series_cols

strip_style <- if (isTRUE(config$strip_colour)) "box" else config$strip_style
if (identical(strip_style, "box")) {
    strip <- ggh4x::strip_themed(
        background_x = lapply(series_cols_hex, function(col) {
            ggplot2::element_rect(
                fill = grDevices::adjustcolor(col, alpha.f = 0.5),
                colour = "grey20", linewidth = 0.4
            )
        }),
        text_x = ggplot2::element_text(face = "bold", colour = "grey10", size = 10)
    )
} else if (identical(strip_style, "text")) {
    strip <- ggh4x::strip_themed(
        background_x = ggplot2::element_blank(),
        text_x = lapply(series_cols_hex, function(col) {
            ggplot2::element_text(colour = col, size = 13)
        })
    )
} else {
    strip <- ggh4x::strip_themed(
        background_x = ggplot2::element_blank(),
        text_x = ggplot2::element_text(face = "bold", size = 10)
    )
}

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
        colour = "#dbdee4",
        linewidth = 0.35
    ) +
    ggplot2::geom_segment(
        data = spokes,
        ggplot2::aes(x = 0, y = 0, xend = px, yend = py),
        colour = "#cfd2d7",
        linewidth = 0.4
    )
if (isTRUE(config$outer_frame)) {
    p <- p + ggplot2::geom_path(
        data = rings[rings$ring == config$n_rings, , drop = FALSE],
        ggplot2::aes(x = px, y = py), colour = "grey15", linewidth = 0.6
    )
}
if (identical(fill_mode, "baseline")) {
    p <- p + ggplot2::geom_polygon(
        data = band_pieces,
        ggplot2::aes(x = px, y = py, group = id, fill = series),
        alpha = config$fill_alpha, colour = NA
    )
    if (isTRUE(config$baseline_line)) {
        p <- p + ggplot2::geom_path(
            data = base_path, ggplot2::aes(x = px, y = py),
            colour = "grey55", linewidth = 0.35
        )
    }
    p <- p + ggplot2::geom_path(
        data = long_closed,
        ggplot2::aes(x = px, y = py, colour = series, group = series),
        linewidth = 0.8, linejoin = "round"
    )
} else {
    p <- p + ggplot2::geom_polygon(
        data = long_closed,
        ggplot2::aes(x = px, y = py, colour = series, fill = series, group = series),
        alpha = config$fill_alpha,
        linewidth = 0.8
    )
}
p <- p +
    (if (identical(config$point_style, "hollow")) {
        ggplot2::geom_point(
            data = long,
            ggplot2::aes(x = px, y = py, colour = series),
            shape = 21, fill = "white", size = 2.2, stroke = 0.9
        )
    } else {
        ggplot2::geom_point(
            data = long,
            ggplot2::aes(x = px, y = py, colour = series),
            size = 2
        )
    }) +
    ggplot2::geom_text(
        data = axis_lab,
        ggplot2::aes(x = px, y = py, label = label, hjust = hjust, vjust = vjust),
        size = config$axis_label_size,
        lineheight = 0.9
    ) +
    ggplot2::scale_colour_manual(values = series_cols_hex, name = NULL) +
    ggplot2::scale_fill_manual(values = series_cols_hex, name = NULL)

if (!is.null(band)) {
    p <- p + ggplot2::geom_text(
        data = band_lab,
        ggplot2::aes(x = px, y = py, label = label, angle = angle),
        size = 2.6,
        colour = "grey30",
        fontface = "italic"
    )
}

if (!is.null(val_lab)) {
    p <- p + ggplot2::geom_text(
        data = val_lab,
        ggplot2::aes(x = px, y = py, label = label),
        size = config$value_label_size,
        fontface = "bold",
        colour = "grey25"
    )
}

lim <- lab_r + 0.05
n_panel <- length(series_cols)
ncol_f <- max(1L, min(as.integer(config$ncol), n_panel))
nrow_f <- ceiling(n_panel / ncol_f)
p <- p +
    ggh4x::facet_wrap2(~series, ncol = ncol_f, strip = strip) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim + 0.12),
        expand = FALSE, clip = "off") +
    ggplot2::labs(title = config$labels$title, x = NULL, y = NULL) +
    theme_prism() +
    ggplot2::theme(
        plot.margin = ggplot2::margin(16, 70, 16, 70),
        plot.background = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, margin = ggplot2::margin(b = 16)),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        panel.spacing.x = ggplot2::unit(42, "mm"),
        panel.spacing.y = ggplot2::unit(16, "mm"),
        strip.text = ggplot2::element_text(margin = ggplot2::margin(4, 4, 4, 4)),
        strip.clip = "off",
        # Series-coloured titles already identify each panel.
        legend.position = if (identical(strip_style, "text")) "none" else "bottom",
        legend.direction = "horizontal"
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 4.6 * ncol_f + 1.6, height = 4.0 * nrow_f + 1.2)
