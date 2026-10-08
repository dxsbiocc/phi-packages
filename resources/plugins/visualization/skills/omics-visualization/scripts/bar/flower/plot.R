#!/usr/bin/env Rscript

# Template-ID: bar-flower
#
# Purpose:
#   Draw a multi-set flower (petal) diagram: a centre core shared by
#   all sets, and one petal per set with supplied total and unique
#   counts. Use when five or more sets make classic Venn unreadable.
#
# Inputs:
#   One row per set. Default example:
#     - set: set / group name
#     - total: items in that set (any overlap)
#     - unique: items exclusive to that set
#   Core count is CONFIG (shared by all sets).
#
# Output:
#   A PDF, PNG, or SVG flower plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggforce
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, core_n, palette, petal geometry,
#   labels). Extra sets are extra rows, not a new id. Do not compute
#   set membership from raw samples here — supply total / unique /
#   core. Two–four overlapping circles stay bar-venn; combination
#   matrices stay bar-upset. Petals are annular sectors (radial
#   sides + outer arc), not ellipses.
#
# Scientific assumptions:
#   total, unique, and core_n are supplied display counts. This
#   script does not enumerate items, test enrichment, or estimate
#   overlaps.

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
        set = "set",
        total = "total",
        unique = "unique"
    ),
    core_n = 936L,
    set_order = NULL,
    # Reference flower colours (CK → Mix-H, clockwise from top).
    # Swap for a named id like "Qualitative.Pastel" if preferred.
    palettes = list(
        petal = c(
            "#BDBDBD", "#F8CDAF", "#F09C87", "#EF645E", "#A3DEC8",
            "#6ECDD6", "#76A6ED", "#ADB8F8", "#C7A5E2", "#F9C2DF"
        )
    ),
    geometry = list(
        core_r = 0.36,
        # Inner / outer radius of each annular petal (gap to core + tip).
        petal_inner = 0.44,
        petal_outer = 1.08,
        # Fraction of each angular slot left empty between petals.
        petal_gap = 0.16,
        # Outer-corner fillet radius (圆角).
        corner_r = 0.15,
        label_r = 1.26,
        n_arc = 48L,
        n_corner = 16L
    ),
    labels = list(
        title = "",
        core = "Core"
    ),
    size = list(
        width = 8.2,
        height = 8.2
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggforce"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

set_col <- config$columns$set
tot_col <- config$columns$total
uni_col <- config$columns$unique

df[[set_col]] <- as.character(df[[set_col]])
df[[tot_col]] <- as.numeric(df[[tot_col]])
df[[uni_col]] <- as.numeric(df[[uni_col]])

ok <- nzchar(df[[set_col]]) &
    is.finite(df[[tot_col]]) & df[[tot_col]] >= 0 &
    is.finite(df[[uni_col]]) & df[[uni_col]] >= 0
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with invalid set/total/unique.")
    df <- df[ok, , drop = FALSE]
}
if (nrow(df) < 3L) {
    stop("Flower plot needs at least 3 sets.", call. = FALSE)
}
if (anyDuplicated(df[[set_col]])) {
    stop("set names must be unique.", call. = FALSE)
}

core_n <- as.integer(config$core_n)
if (!is.finite(core_n) || core_n < 0) {
    stop("config$core_n must be a non-negative integer.", call. = FALSE)
}

ord <- config$set_order
if (is.null(ord) || !length(ord)) {
    ord <- df[[set_col]]
} else {
    miss <- setdiff(ord, df[[set_col]])
    if (length(miss)) {
        stop(
            "set_order missing from table: ",
            paste(miss, collapse = ", "),
            call. = FALSE
        )
    }
    df <- df[match(ord, df[[set_col]]), , drop = FALSE]
}

n <- nrow(df)
df$idx <- seq_len(n)
# Clockwise from top (CK at 12 o'clock in the default example).
df$theta <- pi / 2 - (df$idx - 1L) * (2 * pi / n)
petal_pal <- config$palettes$petal
if (length(petal_pal) == 1L && !startsWith(petal_pal, "#")) {
    df$fill <- palette_colors(petal_pal, n = n)
} else {
    df$fill <- expand_palette(petal_pal, n)
}
df$lab_total <- as.character(as.integer(round(df[[tot_col]])))
df$lab_unique <- paste0("(", as.integer(round(df[[uni_col]])), ")")
df$tx <- ((config$geometry$petal_inner + config$geometry$petal_outer) / 2) *
    cos(df$theta)
df$ty <- ((config$geometry$petal_inner + config$geometry$petal_outer) / 2) *
    sin(df$theta)
df$lx <- config$geometry$label_r * cos(df$theta)
df$ly <- config$geometry$label_r * sin(df$theta)

# Annular-sector petal with filleted outer corners (圆角外侧).
# Radial sides + concentric outer arc; not an ellipse.
make_petal <- function(id, theta, fill, r_in, r_out, half_ang,
                       corner_r, n_arc, n_corner) {
    a0 <- theta - half_ang
    a1 <- theta + half_ang
    # Keep fillet inside the petal band and angular half-width.
    rr <- min(
        corner_r,
        0.48 * (r_out - r_in),
        0.85 * half_ang * (r_out - corner_r)
    )
    rr <- max(rr, 0)

    arc_pts <- function(cx, cy, rad, from, to, n) {
        ang <- seq(from, to, length.out = n)
        cbind(cx + rad * cos(ang), cy + rad * sin(ang))
    }

    if (rr < 1e-4) {
        ao <- seq(a0, a1, length.out = n_arc)
        ai <- seq(a1, a0, length.out = n_arc)
        xy <- rbind(
            cbind(r_out * cos(ao), r_out * sin(ao)),
            cbind(r_in * cos(ai), r_in * sin(ai))
        )
    } else {
        # Angular inset so the fillet is tangent to the outer circle
        # and to each radial side.
        alpha <- asin(rr / (r_out - rr))
        # Left outer fillet centre (toward petal interior from a0).
        phi0 <- a0 + alpha
        c0x <- (r_out - rr) * cos(phi0)
        c0y <- (r_out - rr) * sin(phi0)
        # Right outer fillet centre (toward interior from a1).
        phi1 <- a1 - alpha
        c1x <- (r_out - rr) * cos(phi1)
        c1y <- (r_out - rr) * sin(phi1)

        # Tangent on left radial side.
        t0_r <- (r_out - rr) * cos(alpha)
        # Fillet sweep: from radial tangent to outer-circle tangent.
        # Vectors from centre to tangent points.
        v0_rad <- c(t0_r * cos(a0) - c0x, t0_r * sin(a0) - c0y)
        v0_out <- c(r_out * cos(phi0) - c0x, r_out * sin(phi0) - c0y)
        ang0_from <- atan2(v0_rad[[2]], v0_rad[[1]])
        ang0_to <- atan2(v0_out[[2]], v0_out[[1]])
        # Unwrap short CCW/CW arc staying outside.
        while (ang0_to - ang0_from > pi) ang0_to <- ang0_to - 2 * pi
        while (ang0_to - ang0_from < -pi) ang0_to <- ang0_to + 2 * pi

        t1_r <- (r_out - rr) * cos(alpha)
        v1_out <- c(r_out * cos(phi1) - c1x, r_out * sin(phi1) - c1y)
        v1_rad <- c(t1_r * cos(a1) - c1x, t1_r * sin(a1) - c1y)
        ang1_from <- atan2(v1_out[[2]], v1_out[[1]])
        ang1_to <- atan2(v1_rad[[2]], v1_rad[[1]])
        while (ang1_to - ang1_from > pi) ang1_to <- ang1_to - 2 * pi
        while (ang1_to - ang1_from < -pi) ang1_to <- ang1_to + 2 * pi

        # Left radial (inner -> fillet).
        left_rad <- cbind(
            seq(r_in, t0_r, length.out = 8) * cos(a0),
            seq(r_in, t0_r, length.out = 8) * sin(a0)
        )
        left_fillet <- arc_pts(c0x, c0y, rr, ang0_from, ang0_to, n_corner)
        outer_arc <- cbind(
            r_out * cos(seq(phi0, phi1, length.out = n_arc)),
            r_out * sin(seq(phi0, phi1, length.out = n_arc))
        )
        right_fillet <- arc_pts(c1x, c1y, rr, ang1_from, ang1_to, n_corner)
        right_rad <- cbind(
            seq(t1_r, r_in, length.out = 8) * cos(a1),
            seq(t1_r, r_in, length.out = 8) * sin(a1)
        )
        inner_arc <- cbind(
            r_in * cos(seq(a1, a0, length.out = n_arc)),
            r_in * sin(seq(a1, a0, length.out = n_arc))
        )
        xy <- rbind(
            left_rad, left_fillet, outer_arc,
            right_fillet, right_rad, inner_arc
        )
    }

    data.frame(
        id = id,
        fill = fill,
        x = xy[, 1],
        y = xy[, 2],
        stringsAsFactors = FALSE
    )
}

slot <- 2 * pi / n
half_ang <- slot * (1 - config$geometry$petal_gap) / 2
petal_df <- do.call(
    rbind,
    lapply(seq_len(n), function(i) {
        make_petal(
            id = i,
            theta = df$theta[[i]],
            fill = df$fill[[i]],
            r_in = config$geometry$petal_inner,
            r_out = config$geometry$petal_outer,
            half_ang = half_ang,
            corner_r = config$geometry$corner_r,
            n_arc = as.integer(config$geometry$n_arc),
            n_corner = as.integer(config$geometry$n_corner)
        )
    })
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
core_fill <- "#ffffff"

p <- ggplot() +
    geom_polygon(
        data = petal_df,
        aes(x = x, y = y, group = id, fill = I(fill)),
        colour = NA
    ) +
    ggforce::geom_circle(
        aes(
            x0 = 0, y0 = 0,
            r = config$geometry$core_r
        ),
        fill = core_fill,
        colour = NA,
        inherit.aes = FALSE
    ) +
    geom_text(
        data = df,
        aes(x = tx, y = ty, label = lab_total),
        size = 3.4,
        fontface = "bold",
        colour = "grey15",
        vjust = 0.1
    ) +
    geom_text(
        data = df,
        aes(x = tx, y = ty, label = lab_unique),
        size = 2.8,
        colour = "grey25",
        vjust = 1.4
    ) +
    annotate(
        "text",
        x = 0, y = 0.06,
        label = config$labels$core,
        size = 4.2,
        fontface = "bold",
        colour = "grey10"
    ) +
    annotate(
        "text",
        x = 0, y = -0.1,
        label = as.character(core_n),
        size = 4.6,
        fontface = "bold",
        colour = "grey10"
    ) +
    geom_text(
        data = df,
        aes(x = lx, y = ly, label = .data[[set_col]]),
        size = 3.6,
        fontface = "bold",
        colour = "grey10"
    ) +
    coord_fixed(clip = "off") +
    labs(title = config$labels$title) +
    theme_void() +
    theme(
        plot.background = element_blank(),
        plot.title = element_text(
            face = "bold", size = 12, hjust = 0.5,
            margin = margin(b = 6)
        ),
        plot.margin = margin(12, 18, 12, 18)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p, io$output,
    width = config$size$width,
    height = config$size$height
)
