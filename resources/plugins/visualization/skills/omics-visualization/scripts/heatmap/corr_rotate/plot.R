#!/usr/bin/env Rscript

# Template-ID: heatmap-corr-rotate
#
# Purpose:
#   Draw a lower-triangle correlation heatmap rotated 45° so the
#   matrix reads as an upright triangle. Optional variable groups
#   insert white gaps; optional annotation tracks sit along the
#   left hypotenuse; non-significant cells can show a "+" mark.
#
# Inputs:
#   1) Long correlation table (one row per pair):
#        x, y, rho[, p]
#   2) Optional annotation table (one row per variable):
#        var, [Group], <track columns...>
#
# Output:
#   A PDF, PNG, or SVG rotated triangular correlation heatmap.
#
# Dependencies:
#   ggplot2, readr, ggprism, scales
#
# Example:
#   Rscript plot.R example.tsv annot.tsv output.pdf
#   Rscript plot.R example.tsv NA output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, annotation tracks, group gap, p cutoff,
#   palette, and labels. Rotation and triangle side stay in CONFIG —
#   not a new id. Do not compute Spearman here. A flat lower triangle
#   with star cells is heatmap-shape; braced flat groups are
#   heatmap-corr-grouped.
#
# Scientific assumptions:
#   rho (and optional p) are supplied display values. This script
#   does not adjust p-values or impute missing pairs.

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

io <- parse_io_args(
    n_input = 2L,
    input_names = c("corr", "annot")
)

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        x = "x",
        y = "y",
        value = "rho",
        p = "p",
        var = "var",
        group = "Group"
    ),
    # Annotation track columns from the second table (after var/Group).
    annot_tracks = c("Detection", "Species"),
    # Gap inserted between adjacent variable groups (index units).
    group_gap = 0.35,
    # p > ns_p → grey tile with "+" (set NULL to disable).
    ns_p = 0.05,
    var_order = NULL,
    palettes = list(
        # Spectral (blue/purple ↔ red) matches the reference figure.
        fill = "Diverging.Spectral / 11",
        # Fixed legend colours from the reference (Pos/Neg, species).
        Detection = c(Pos = "#5B8FD9", Neg = "#E87A8A"),
        Species = c(
            `spec-A` = "#E89A3C",
            `spec-B` = "#5BA86B",
            `spec-C` = "#4A7AB5"
        )
    ),
    fill_limits = c(-1, 1),
    geometry = list(
        tile_s = 0.40,
        corner_r = 0.08,
        n_corner = 8L,
        border = FALSE,
        group_outline = TRUE,
        outline_pad = 0.12
    ),
    labels = list(
        title = "",
        fill = "Correlation",
        x_size = 9
    ),
    size = list(
        width = 9.6,
        height = 7.2
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "scales", "sf"))

# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
# Axis-aligned rounded square → rotate 45° → rounded diamond in (u, v).
rounded_diamond <- function(u, v, s = 0.42, rr = 0.11, n_arc = 8L) {
    rr <- min(rr, 0.45 * s)
    half <- s / sqrt(2)
    # Rounded square at origin, then rotate 45° → diamond.
    arc <- function(cx, cy, from, to, n) {
        ang <- seq(from, to, length.out = n)
        cbind(cx + rr * cos(ang), cy + rr * sin(ang))
    }
    sq <- rbind(
        arc(half - rr, half - rr, 0, pi / 2, n_arc),
        arc(-half + rr, half - rr, pi / 2, pi, n_arc),
        arc(-half + rr, -half + rr, pi, 3 * pi / 2, n_arc),
        arc(half - rr, -half + rr, 3 * pi / 2, 2 * pi, n_arc)
    )
    rot <- pi / 4
    xr <- sq[, 1] * cos(rot) - sq[, 2] * sin(rot)
    yr <- sq[, 1] * sin(rot) + sq[, 2] * cos(rot)
    data.frame(x = u + xr, y = v + yr)
}

# Outer contour of a tile set: buffer each polygon by `pad`, then union.
# Uses sf so the outline is a true Euclidean offset of the tile shapes
# (zigzag base, rounded corners, uniform clearance).
outline_from_tiles <- function(tile_df, pad, fid) {
    ids <- unique(tile_df$id)
    if (!length(ids)) {
        return(NULL)
    }
    geoms <- lapply(ids, function(i) {
        xy <- as.matrix(tile_df[tile_df$id == i, c("x", "y")])
        if (nrow(xy) < 3L) {
            return(NULL)
        }
        if (!isTRUE(all.equal(xy[1, ], xy[nrow(xy), ], tolerance = 1e-9))) {
            xy <- rbind(xy, xy[1, , drop = FALSE])
        }
        sf::st_polygon(list(xy))
    })
    geoms <- Filter(Negate(is.null), geoms)
    if (!length(geoms)) {
        return(NULL)
    }
    g <- sf::st_make_valid(sf::st_sfc(geoms))
    # Buffer then union: merges neighbours when gap < 2*pad, keeps
    # concave zigzags where the notch is deeper than the pad.
    uni <- sf::st_union(sf::st_buffer(g, dist = pad, nQuadSegs = 8L))
    uni <- sf::st_make_valid(uni)
    gtype <- as.character(sf::st_geometry_type(uni, by_geometry = FALSE))
    if (grepl("MULTI", gtype)) {
        parts <- sf::st_cast(uni, "POLYGON")
        areas <- as.numeric(sf::st_area(parts))
        uni <- parts[which.max(areas)]
    }
    coords <- sf::st_coordinates(uni)
    if (is.null(coords) || nrow(coords) < 3L) {
        return(NULL)
    }
    if ("L1" %in% colnames(coords)) {
        coords <- coords[coords[, "L1"] == 1, , drop = FALSE]
    }
    if (nrow(coords) < 3L) {
        return(NULL)
    }
    data.frame(
        id = fid,
        x = as.numeric(coords[, "X"]),
        y = as.numeric(coords[, "Y"]),
        stringsAsFactors = FALSE
    )
}

resolve_track_colors <- function(pal, levels) {
    levels <- as.character(levels)
    if (is.null(pal)) {
        cols <- palette_colors("Qualitative.Safe", n = max(length(levels), 1L))
        return(setNames(cols[seq_along(levels)], levels))
    }
    if (length(pal) == 1L && !startsWith(as.character(pal), "#") &&
        is.null(names(pal))) {
        cols <- palette_colors(pal, n = max(length(levels), 1L))
        return(setNames(cols[seq_along(levels)], levels))
    }
    named <- pal
    if (is.null(names(named)) || !any(nzchar(names(named)))) {
        named <- setNames(as.character(pal), levels[seq_along(pal)])
    }
    out <- character(length(levels))
    names(out) <- levels
    for (lv in levels) {
        if (lv %in% names(named)) {
            out[[lv]] <- unname(named[[lv]])
        } else {
            out[[lv]] <- NA_character_
        }
    }
    miss <- is.na(out)
    if (any(miss)) {
        fill <- palette_colors("Qualitative.Safe", n = sum(miss))
        out[miss] <- fill
    }
    out
}

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
corr <- read_table_auto(io$input[[1]])
annot_path <- io$input[[2]]
has_annot <- !(is.na(annot_path) ||
    identical(toupper(as.character(annot_path)), "NA") ||
    !nzchar(annot_path) ||
    !file.exists(annot_path))

cx <- config$columns$x
cy <- config$columns$y
cv <- config$columns$value
cp <- config$columns$p
gv <- config$columns$var
gg <- config$columns$group

require_columns(corr, list(x = cx, y = cy, value = cv))
corr[[cx]] <- as.character(corr[[cx]])
corr[[cy]] <- as.character(corr[[cy]])
corr[[cv]] <- as.numeric(corr[[cv]])
if (!is.null(cp) && cp %in% names(corr)) {
    corr[[cp]] <- as.numeric(corr[[cp]])
} else {
    cp <- NULL
}

ann <- NULL
if (has_annot) {
    ann <- read_table_auto(annot_path)
    require_columns(ann, list(var = gv))
    ann[[gv]] <- as.character(ann[[gv]])
}

vars <- unique(c(corr[[cx]], corr[[cy]]))
if (!is.null(config$var_order) && length(config$var_order)) {
    vars <- intersect(config$var_order, vars)
} else if (!is.null(ann)) {
    ord <- intersect(ann[[gv]], vars)
    if (length(ord)) {
        vars <- c(ord, setdiff(vars, ord))
    }
}
if (length(vars) < 3L) {
    stop("Need at least 3 variables for a rotated triangle.", call. = FALSE)
}

# Group membership (for gaps + outlines).
grp <- rep("all", length(vars))
names(grp) <- vars
if (!is.null(ann) && !is.null(gg) && gg %in% names(ann)) {
    gmap <- setNames(as.character(ann[[gg]]), ann[[gv]])
    grp <- ifelse(vars %in% names(gmap), unname(gmap[vars]), grp)
    names(grp) <- vars
}

# Positions along the base with gaps between groups.
gap <- as.numeric(config$group_gap)
if (!is.finite(gap) || gap < 0) gap <- 0
pos <- numeric(length(vars))
pos[[1]] <- 1
for (k in seq_along(vars)[-1]) {
    pos[[k]] <- pos[[k - 1]] + 1 +
        if (!identical(grp[[k]], grp[[k - 1]])) gap else 0
}
names(pos) <- vars

idx <- setNames(seq_along(vars), vars)
corr <- corr[corr[[cx]] %in% vars & corr[[cy]] %in% vars, , drop = FALSE]
corr$i <- unname(idx[corr[[cx]]])
corr$j <- unname(idx[corr[[cy]]])
# Keep lower triangle including diagonal (j >= i).
corr <- corr[corr$j >= corr$i, , drop = FALSE]
ok <- is.finite(corr[[cv]])
if (any(!ok)) {
    message("Dropped ", sum(!ok), " non-finite rho row(s).")
    corr <- corr[ok, , drop = FALSE]
}

corr$pi <- unname(pos[corr[[cx]]])
corr$pj <- unname(pos[corr[[cy]]])
# Rotate 45°: diagonal becomes the base of an upright triangle.
corr$u <- (corr$pi + corr$pj) / 2
corr$v <- (corr$pj - corr$pi) / 2
corr$group_i <- unname(grp[corr[[cx]]])
corr$group_j <- unname(grp[corr[[cy]]])

ns_p <- config$ns_p
corr$ns <- FALSE
if (!is.null(ns_p) && !is.null(cp) && is.finite(ns_p)) {
    corr$ns <- is.finite(corr[[cp]]) &
        corr[[cp]] > ns_p &
        corr$i != corr$j
}

fill_lim <- as.numeric(config$fill_limits)
if (length(fill_lim) != 2L || any(!is.finite(fill_lim))) {
    fill_lim <- c(-1, 1)
}
# Spectral / 11 is red → purple; low rho = red, high rho = blue/purple.
fill_cols <- palette_colors(config$palettes$fill)

lab_pad <- config$geometry$outline_pad
if (is.null(lab_pad) || !is.finite(as.numeric(lab_pad))) {
    lab_pad <- 0.12
}
lab <- data.frame(
    var = vars,
    u = unname(pos),
    v = -(config$geometry$tile_s + as.numeric(lab_pad) + 0.22),
    stringsAsFactors = FALSE
)

# Annotation tracks along the left hypotenuse (outward).
# Positions follow `pos` so group gaps match the heatmap.
ann_tiles <- NULL
ann_guides <- list()
n_var <- length(vars)
apex <- c((pos[[1]] + pos[[n_var]]) / 2, (pos[[n_var]] - pos[[1]]) / 2)
left0 <- c(pos[[1]], 0)
right0 <- c(pos[[n_var]], 0)
edge <- apex - left0
elen <- sqrt(sum(edge^2))
if (elen < 1e-8) elen <- 1
tang <- edge / elen
out_n <- c(-tang[[2]], tang[[1]])
if (out_n[[1]] > 0) out_n <- -out_n
# Right hypotenuse (mirror of left): outward normal points right.
r_edge <- apex - right0
r_elen <- sqrt(sum(r_edge^2))
if (r_elen < 1e-8) r_elen <- 1
r_tang <- r_edge / r_elen
r_out <- c(-r_tang[[2]], r_tang[[1]])
if (r_out[[1]] < 0) r_out <- -r_out
# Scale: one unit of pos → distance along the left edge.
pos_span <- max(pos[[n_var]] - pos[[1]], 1e-8)
pos_to_s <- function(p) (p - pos[[1]]) / pos_span * elen

# Correlation colorbar along the right hypotenuse (like left annots).
cbar_poly <- NULL
cbar_ticks <- NULL
cbar_title <- NULL
{
    pad <- config$geometry$outline_pad
    if (is.null(pad) || !is.finite(as.numeric(pad))) pad <- 0.12
    unit_s <- pos_to_s(pos[[1]] + 1) - pos_to_s(pos[[1]])
    bar_half_w <- unit_s * 0.36
    # Same standoff family as left annotation tracks.
    bar_off <- as.numeric(config$geometry$tile_s) + 0.55
    s0 <- unit_s * 0.08
    s1 <- r_elen - unit_s * 0.08
    n_seg <- 100L
    ss <- seq(s0, s1, length.out = n_seg + 1L)
    # High correlation at apex; low at base (matches tile fill).
    t_mid <- (ss[-length(ss)] + ss[-1]) / 2
    t_frac <- (t_mid - s0) / max(s1 - s0, 1e-8)
    rho_mid <- fill_lim[[1]] + t_frac * (fill_lim[[2]] - fill_lim[[1]])
    pal_fun <- scales::col_numeric(
        palette = fill_cols,
        domain = fill_lim,
        na.color = "#E8E8E8"
    )
    cols <- pal_fun(rho_mid)
    origin <- right0 + bar_off * r_out
    rows <- lapply(seq_len(n_seg), function(i) {
        s_a <- ss[[i]]
        s_b <- ss[[i + 1L]]
        p1 <- origin + s_a * r_tang
        p2 <- origin + s_b * r_tang
        corners <- rbind(
            p1 + bar_half_w * r_out,
            p2 + bar_half_w * r_out,
            p2 - bar_half_w * r_out,
            p1 - bar_half_w * r_out
        )
        data.frame(
            id = i,
            fill = cols[[i]],
            x = corners[, 1],
            y = corners[, 2],
            stringsAsFactors = FALSE
        )
    })
    cbar_poly <- do.call(rbind, rows)
    brks <- seq(fill_lim[[1]], fill_lim[[2]], by = 0.2)
    tick_s <- s0 + (brks - fill_lim[[1]]) /
        max(fill_lim[[2]] - fill_lim[[1]], 1e-8) * (s1 - s0)
    tick_off <- bar_half_w + unit_s * 0.28
    cbar_ticks <- data.frame(
        lab = format(brks, trim = TRUE),
        x = origin[[1]] + tick_s * r_tang[[1]] + tick_off * r_out[[1]],
        y = origin[[2]] + tick_s * r_tang[[2]] + tick_off * r_out[[2]],
        stringsAsFactors = FALSE
    )
    mid_s <- (s0 + s1) / 2
    title_off <- bar_half_w + unit_s * 1.35
    cbar_title <- data.frame(
        lab = config$labels$fill,
        x = origin[[1]] + mid_s * r_tang[[1]] + title_off * r_out[[1]],
        y = origin[[2]] + mid_s * r_tang[[2]] + title_off * r_out[[2]],
        stringsAsFactors = FALSE
    )
}
if (!is.null(ann)) {
    tracks <- intersect(config$annot_tracks, names(ann))
    if (length(tracks)) {
        # Chip span along edge from midpoints in pos-space (keeps
        # within-group chips abutting; between-group gap stays empty).
        s_lo <- numeric(n_var)
        s_hi <- numeric(n_var)
        for (i in seq_len(n_var)) {
            # Within a group: abut at midpoints. At a group break: stop
            # at ±0.5 so the group_gap stays empty on the annot strip.
            same_prev <- i > 1L && identical(grp[[i]], grp[[i - 1]])
            same_next <- i < n_var && identical(grp[[i]], grp[[i + 1]])
            p_lo <- if (same_prev) {
                (pos[[i - 1]] + pos[[i]]) / 2
            } else {
                pos[[i]] - 0.5
            }
            p_hi <- if (same_next) {
                (pos[[i]] + pos[[i + 1]]) / 2
            } else {
                pos[[i]] + 0.5
            }
            s_lo[[i]] <- pos_to_s(p_lo) + 0.004
            s_hi[[i]] <- pos_to_s(p_hi) - 0.004
        }
        rows <- list()
        # Square-ish chips: thickness ≈ within-group half-step on edge.
        unit_s <- pos_to_s(pos[[1]] + 1) - pos_to_s(pos[[1]])
        track_half_w <- unit_s * 0.42
        track_gap <- unit_s * 0.08
        for (ti in seq_along(tracks)) {
            tr <- tracks[[ti]]
            vals <- as.character(ann[[tr]][match(vars, ann[[gv]])])
            levs <- unique(vals[!is.na(vals) & nzchar(vals)])
            col_map <- resolve_track_colors(config$palettes[[tr]], levs)
            off <- config$geometry$tile_s + 0.55 +
                (ti - 1L) * (2 * track_half_w + track_gap)
            u0 <- left0[[1]] + off * out_n[[1]]
            v0 <- left0[[2]] + off * out_n[[2]]
            for (i in seq_len(n_var)) {
                sc <- (s_lo[[i]] + s_hi[[i]]) / 2
                hl <- max((s_hi[[i]] - s_lo[[i]]) / 2, 0.04)
                rows[[length(rows) + 1L]] <- data.frame(
                    track = tr,
                    var = vars[[i]],
                    u = u0 + sc * tang[[1]],
                    v = v0 + sc * tang[[2]],
                    half_len = hl,
                    half_w = track_half_w,
                    level = vals[[i]],
                    fill = unname(col_map[[vals[[i]]]]),
                    stringsAsFactors = FALSE
                )
            }
            ann_guides[[tr]] <- col_map
        }
        ann_tiles <- do.call(rbind, rows)
    }
}

# -----------------------------------------------------------------------------
# TILE GEOMETRY
# -----------------------------------------------------------------------------
gs <- config$geometry
tile_s <- gs$tile_s
corner_r <- gs$corner_r
n_corner <- as.integer(gs$n_corner)

tile_list <- lapply(seq_len(nrow(corr)), function(k) {
    d <- rounded_diamond(
        corr$u[[k]], corr$v[[k]],
        s = tile_s, rr = corner_r, n_arc = n_corner
    )
    d$id <- k
    d$rho <- corr[[cv]][[k]]
    d$ns <- corr$ns[[k]]
    d$group_i <- corr$group_i[[k]]
    d$group_j <- corr$group_j[[k]]
    d
})
tiles <- do.call(rbind, tile_list)
tiles_sig <- tiles[!tiles$ns, , drop = FALSE]
tiles_ns <- tiles[tiles$ns, , drop = FALSE]

# Block outlines via sf: buffer each tile polygon, union → Euclidean contour.
outline_df <- NULL
if (isTRUE(gs$group_outline) && length(unique(grp)) > 1L) {
    ug <- unique(grp)
    pad <- as.numeric(gs$outline_pad)
    if (!is.finite(pad) || pad <= 0) {
        pad <- 0.12
    }
    frames <- list()
    fid <- 0L

    # Within-group triangles.
    for (gname in ug) {
        sub <- tiles[tiles$group_i == gname & tiles$group_j == gname, , drop = FALSE]
        if (!nrow(sub)) next
        fid <- fid + 1L
        frames[[fid]] <- outline_from_tiles(sub, pad, fid)
    }

    # Between-group rectangles.
    if (length(ug) >= 2L) {
        for (a in seq_len(length(ug) - 1L)) {
            for (b in seq(a + 1L, length(ug))) {
                ga <- ug[[a]]
                gb <- ug[[b]]
                sub <- tiles[
                    (tiles$group_i == ga & tiles$group_j == gb) |
                        (tiles$group_i == gb & tiles$group_j == ga),
                    ,
                    drop = FALSE
                ]
                if (!nrow(sub)) next
                fid <- fid + 1L
                frames[[fid]] <- outline_from_tiles(sub, pad, fid)
            }
        }
    }
    frames <- Filter(Negate(is.null), frames)
    if (length(frames)) {
        outline_df <- do.call(rbind, frames)
    }
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
border_col <- if (isTRUE(gs$border)) "grey90" else NA
border_lw <- if (isTRUE(gs$border)) 0.25 else 0

p <- ggplot()

if (nrow(tiles_sig)) {
    p <- p +
        geom_polygon(
            data = tiles_sig,
            aes(x = x, y = y, group = id, fill = rho),
            colour = border_col,
            linewidth = border_lw
        )
}

if (nrow(tiles_ns)) {
    p <- p +
        geom_polygon(
            data = tiles_ns,
            aes(x = x, y = y, group = id),
            fill = "#E8E8E8",
            colour = border_col,
            linewidth = border_lw,
            inherit.aes = FALSE
        ) +
        geom_text(
            data = corr[corr$ns, , drop = FALSE],
            aes(x = u, y = v),
            label = "+",
            size = 3.2,
            colour = "grey20",
            inherit.aes = FALSE
        )
}

# Frames: sf buffer∪union contour of each block (uniform pad).
if (!is.null(outline_df) && nrow(outline_df)) {
    p <- p +
        geom_polygon(
            data = outline_df,
            aes(x = x, y = y, group = id),
            fill = NA,
            colour = "grey15",
            linewidth = 0.55,
            linejoin = "round",
            inherit.aes = FALSE
        )
}

p <- p +
    scale_fill_gradientn(
        name = config$labels$fill,
        colours = fill_cols,
        limits = fill_lim,
        breaks = seq(fill_lim[[1]], fill_lim[[2]], by = 0.2),
        oob = scales::squish,
        guide = "none"
    )

if (!is.null(ann_tiles) && nrow(ann_tiles)) {
    # Rectangles: long side along the hypotenuse, short side along
    # the outward normal; neighbours share edges (tiny gap_pad).
    ann_poly <- do.call(rbind, lapply(seq_len(nrow(ann_tiles)), function(k) {
        u <- ann_tiles$u[[k]]
        v <- ann_tiles$v[[k]]
        hl <- ann_tiles$half_len[[k]]
        hw <- ann_tiles$half_w[[k]]
        corners <- rbind(
            c(u, v) + hl * tang + hw * out_n,
            c(u, v) - hl * tang + hw * out_n,
            c(u, v) - hl * tang - hw * out_n,
            c(u, v) + hl * tang - hw * out_n
        )
        data.frame(
            id = k,
            fill = ann_tiles$fill[[k]],
            x = corners[, 1],
            y = corners[, 2],
            stringsAsFactors = FALSE
        )
    }))
    p <- p +
        geom_polygon(
            data = ann_poly,
            aes(x = x, y = y, group = id),
            fill = ann_poly$fill,
            colour = NA,
            inherit.aes = FALSE
        )
}

# Slanted correlation colorbar on the right hypotenuse.
if (!is.null(cbar_poly) && nrow(cbar_poly)) {
    # Keep labels upright-along-bar; flip 180° if needed for reading.
    edge_ang <- atan2(r_tang[[2]], r_tang[[1]]) * 180 / pi
    if (edge_ang > 90) edge_ang <- edge_ang - 180
    if (edge_ang < -90) edge_ang <- edge_ang + 180
    p <- p +
        geom_polygon(
            data = cbar_poly,
            aes(x = x, y = y, group = id),
            fill = cbar_poly$fill,
            colour = NA,
            inherit.aes = FALSE
        ) +
        geom_text(
            data = cbar_ticks,
            aes(x = x, y = y, label = lab),
            angle = edge_ang,
            hjust = 0,
            vjust = 0.5,
            size = 2.6,
            colour = "grey10",
            inherit.aes = FALSE
        ) +
        geom_text(
            data = cbar_title,
            aes(x = x, y = y, label = lab),
            angle = edge_ang,
            hjust = 0.5,
            vjust = 0.5,
            fontface = "bold",
            size = 3,
            colour = "grey10",
            inherit.aes = FALSE
        )
}

p <- p +
    geom_text(
        data = lab,
        aes(x = u, y = v, label = var),
        angle = 45,
        hjust = 1,
        vjust = 0.5,
        size = config$labels$x_size / 2.8,
        colour = "grey10",
        inherit.aes = FALSE
    ) +
    coord_fixed(clip = "off") +
    labs(title = config$labels$title) +
    theme_void() +
    theme(
        plot.background = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(
            face = "bold", size = 12, hjust = 0.5
        ),
        plot.margin = margin(24, 56, 30, 36),
        legend.position = "none"
    )

# Manual discrete legends: mid-right of the triangle, outside the
# slanted correlation bar.
if (length(ann_guides) >= 1L) {
    n_title <- length(ann_guides)
    n_keys <- sum(vapply(ann_guides, length, integer(1)))
    # Match the yi steps used below (0.38 per row + 0.35 between groups).
    stack_h <- (n_title + n_keys) * 0.38 + n_title * 0.35
    unit_s <- pos_to_s(pos[[1]] + 1) - pos_to_s(pos[[1]])
    # Past colorbar body, tick labels, and "Correlation" title.
    leg_off <- as.numeric(config$geometry$tile_s) + 0.55 +
        unit_s * 0.36 + unit_s * 4.0
    mid_edge <- (apex + right0) / 2
    # Center on mid-edge, then nudge up so the block reads as mid-right.
    x0 <- mid_edge[[1]] + leg_off * r_out[[1]]
    y0 <- mid_edge[[2]] + stack_h / 2 + 1.1
    leg_rows <- list()
    yi <- 0
    for (tr in names(ann_guides)) {
        cmap <- ann_guides[[tr]]
        leg_rows[[length(leg_rows) + 1L]] <- data.frame(
            kind = "title",
            lab = tr,
            x = x0,
            y = y0 - yi * 0.38,
            fill = NA_character_,
            stringsAsFactors = FALSE
        )
        yi <- yi + 1
        for (lv in names(cmap)) {
            leg_rows[[length(leg_rows) + 1L]] <- data.frame(
                kind = "key",
                lab = lv,
                x = x0,
                y = y0 - yi * 0.38,
                fill = unname(cmap[[lv]]),
                stringsAsFactors = FALSE
            )
            yi <- yi + 1
        }
        yi <- yi + 0.35
    }
    leg <- do.call(rbind, leg_rows)
    keys <- leg[leg$kind == "key", , drop = FALSE]
    tits <- leg[leg$kind == "title", , drop = FALSE]
    if (nrow(keys)) {
        key_poly <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
            u <- keys$x[[k]] - 0.35
            v <- keys$y[[k]]
            w <- 0.14
            h <- 0.14
            data.frame(
                id = k,
                fill = keys$fill[[k]],
                x = c(u - w, u + w, u + w, u - w),
                y = c(v - h, v - h, v + h, v + h),
                stringsAsFactors = FALSE
            )
        }))
        p <- p +
            geom_polygon(
                data = key_poly,
                aes(x = x, y = y, group = id),
                fill = key_poly$fill,
                colour = "grey30",
                linewidth = 0.3,
                inherit.aes = FALSE
            ) +
            geom_text(
                data = keys,
                aes(x = x, y = y, label = lab),
                hjust = 0,
                size = 2.8,
                colour = "grey10",
                inherit.aes = FALSE
            )
    }
    if (nrow(tits)) {
        p <- p +
            geom_text(
                data = tits,
                aes(x = x - 0.35, y = y, label = lab),
                hjust = 0,
                fontface = "bold",
                size = 3,
                colour = "grey10",
                inherit.aes = FALSE
            )
    }
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p, io$output,
    width = config$size$width,
    height = config$size$height
)
