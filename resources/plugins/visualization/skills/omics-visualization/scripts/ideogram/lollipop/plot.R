#!/usr/bin/env Rscript

# Template-ID: ideogram-lollipop
#
# Purpose:
#   Draw a protein mutation lollipop: stems at amino-acid positions,
#   pie or solid-circle heads by mutation class, and a domain bar.
#
# Inputs:
#   One row per mutation or site × class. Default example:
#     - aa: amino-acid position
#     - class: mutation class (already mapped)
#     - n: supplied count (optional; default 1)
#     - label: optional hotspot label
#   Optional sidecar beside the input:
#     - domains.tsv: start, end, domain; optional group (feature colour)
#
# Layouts:
#   layout = "linear" draws the protein axis with a domain bar below.
#   layout = "circular" bends a circular sequence (mtDNA, plasmid, viral
#   genome) into a feature ring, with stems pointing outward (or inward)
#   to solid heads, position labels, and dashed score rings
#   (example_mito.tsv + domains_mito.tsv). Heads are solid circles here.
#
# Output:
#   A PDF, PNG, or SVG protein lollipop figure.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggrepel, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, protein length, pie vs
#   circle, pie radius, label cut). Themes, zoom, and MAF / cBioPortal / Pfam
#   helpers are upstream or CONFIG, not a second id. Transcript exon
#   models are ideogram-gene. Signed or ranked score lollipops are
#   scatter-lollipop-*. Site × amino-acid ΔΔG tiles are
#   heatmap-mutation-energy.
#
# Scientific assumptions:
#   Position, class, and count are supplied. This script does not parse
#   MAF, fetch cBioPortal, look up Pfam, or map Variant_Classification.
#   Summing n at a site is a display aggregate of supplied rows.
#   Domain intervals are a sidecar, not a domain prediction.

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
        pos = "aa",
        class = "class",
        count = "n",
        label = "label"
    ),
    sidecars = list(
        domains = "domains.tsv"
    ),
    domain = list(
        start = "start",
        end = "end",
        name = "domain",
        group = "group"         # optional; colours features by group
    ),
    layout = "linear",          # "linear" or "circular"
    label_prefix = "p.",        # auto labels: prefix + position
    label_big_mark = "",        # e.g. "," gives 12,417
    protein_length = 393,
    pop = "pie",
    pop_r_mm = c(1.0, 3.0),
    label_min = 8,
    fig = list(
        width = 8.8,
        height = 4.2
    ),
    circular = list(
        start_degree = 0,       # sequence position 1 at 12 o'clock
        clockwise = TRUE,
        stem_direction = "outward",  # or "inward"
        stem_length = 1.3,      # longest stem, in ring radii
        ring_width = 0.1,       # feature ring thickness, in ring radii
        rings = NULL,           # dashed score rings; NULL = pretty breaks
        head_size = 3,
        label_size = 3.4,
        legend_ncol = 4,
        width = 7.6,
        height = 7.3
    ),
    labels = list(
        title = NULL,
        x = "Amino acid",
        y = "Count",
        fill = "Class"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggrepel", "patchwork"))

# Pies in millimetres so they stay circular. A full 2π slice is a
# filled circle with no radial stroke (that white radius is not a slice).
GeomPieMm <- ggplot2::ggproto(
    "GeomPieMm",
    ggplot2::Geom,
    required_aes = c("x", "y", "from", "to"),
    default_aes = ggplot2::aes(
        fill = "grey50",
        colour = "white",
        linewidth = 0.15,
        alpha = 1,
        r_mm = 2.8,
        outline = "grey25"
    ),
    extra_params = c("na.rm"),
    draw_key = ggplot2::draw_key_polygon,
    draw_panel = function(data, panel_params, coord) {
        data <- ggplot2::remove_missing(
            data,
            na.rm = TRUE,
            vars = c("x", "y", "from", "to", "r_mm"),
            name = "geom_pie_mm"
        )
        if (!nrow(data)) {
            return(ggplot2::zeroGrob())
        }
        coords <- coord$transform(data, panel_params)
        n <- nrow(coords)
        wedges <- vector("list", n)
        for (i in seq_len(n)) {
            row <- coords[i, ]
            from <- row$from
            to <- row$to
            if (to < from) {
                to <- to + 2 * pi
            }
            span <- to - from
            cx <- grid::unit(row$x, "native")
            cy <- grid::unit(row$y, "native")
            rr <- grid::unit(row$r_mm, "mm")
            if (span >= (2 * pi - 1e-6)) {
                wedges[[i]] <- grid::circleGrob(
                    x = cx,
                    y = cy,
                    r = rr,
                    gp = grid::gpar(fill = row$fill, col = NA, alpha = row$alpha)
                )
            } else {
                n_arc <- max(8L, as.integer(36 * span / (2 * pi)))
                th <- seq(from, to, length.out = n_arc)
                wedges[[i]] <- grid::polygonGrob(
                    x = grid::unit.c(cx, cx + rr * sin(th), cx),
                    y = grid::unit.c(cy, cy + rr * cos(th), cy),
                    gp = grid::gpar(
                        fill = row$fill,
                        col = row$colour,
                        lwd = row$linewidth * ggplot2::.pt,
                        alpha = row$alpha
                    )
                )
            }
        }
        sites <- unique(coords[, c("x", "y", "r_mm", "outline"), drop = FALSE])
        rings <- vector("list", nrow(sites))
        for (i in seq_len(nrow(sites))) {
            rings[[i]] <- grid::circleGrob(
                x = grid::unit(sites$x[[i]], "native"),
                y = grid::unit(sites$y[[i]], "native"),
                r = grid::unit(sites$r_mm[[i]], "mm"),
                gp = grid::gpar(
                    fill = NA,
                    col = sites$outline[[i]],
                    lwd = 0.35 * ggplot2::.pt
                )
            )
        }
        grid::gTree(children = do.call(grid::gList, c(wedges, rings)))
    }
)

geom_pie_mm <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       ..., na.rm = FALSE, inherit.aes = TRUE) {
    ggplot2::layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomPieMm,
        position = position,
        inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...)
    )
}

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
need <- config$columns[c("pos", "class")]
require_columns(df, need)

pos_col <- config$columns$pos
cl_col <- config$columns$class
n_col <- config$columns$count
lab_col <- config$columns$label

df[[pos_col]] <- as.numeric(df[[pos_col]])
if (is.null(n_col) || !nzchar(n_col) || !(n_col %in% names(df))) {
    df$n <- 1
    n_col <- "n"
} else {
    df[[n_col]] <- as.numeric(df[[n_col]])
}
ok <- is.finite(df[[pos_col]]) & df[[pos_col]] > 0 &
    is.finite(df[[n_col]]) & df[[n_col]] > 0 &
    !is.na(df[[cl_col]]) & nzchar(as.character(df[[cl_col]]))
if (any(!ok)) {
    message("Dropped ", sum(!ok), " mutation row(s) with invalid position or count.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No mutation rows left to plot.", call. = FALSE)
}

preferred <- c("Missense", "Truncating", "Inframe", "Other")
present <- unique(as.character(df[[cl_col]]))
class_levels <- c(preferred[preferred %in% present], setdiff(present, preferred))
df[[cl_col]] <- factor(as.character(df[[cl_col]]), levels = class_levels)

agg <- stats::aggregate(
    df[[n_col]],
    by = list(pos = df[[pos_col]], class = df[[cl_col]]),
    FUN = sum
)
names(agg)[3] <- "n"
stems <- stats::aggregate(n ~ pos, data = agg, FUN = sum)
names(stems)[2] <- "total"

lab_by_pos <- NULL
if (!is.null(lab_col) && lab_col %in% names(df)) {
    raw_lab <- trimws(as.character(df[[lab_col]]))
    keep_lab <- !is.na(raw_lab) & nzchar(raw_lab)
    if (any(keep_lab)) {
        lab_by_pos <- stats::aggregate(
            raw_lab[keep_lab],
            by = list(pos = df[[pos_col]][keep_lab]),
            FUN = function(z) z[[1]]
        )
        names(lab_by_pos)[2] <- "label"
    }
}
stems$label <- NA_character_
if (!is.null(lab_by_pos)) {
    stems$label <- lab_by_pos$label[match(stems$pos, lab_by_pos$pos)]
}
need_lab <- is.na(stems$label) & stems$total >= config$label_min
stems$label[need_lab] <- paste0(config$label_prefix,
    format(stems$pos[need_lab], big.mark = config$label_big_mark, trim = TRUE))
stems$label[!(stems$total >= config$label_min)] <- NA_character_

L <- config$protein_length
if (is.null(L) || !is.finite(L) || L < 1) {
    L <- max(stems$pos, na.rm = TRUE)
}

side_path <- function(name) {
    if (is.null(name) || !nzchar(name)) {
        return(NA_character_)
    }
    file.path(dirname(io$input), name)
}

domains <- NULL
dom_path <- side_path(config$sidecars$domains)
if (!is.na(dom_path) && file.exists(dom_path)) {
    domains <- read_table_auto(dom_path)
    require_columns(domains, config$domain[c("start", "end", "name")])
    ds <- config$domain$start
    de <- config$domain$end
    dn <- config$domain$name
    domains[[ds]] <- as.numeric(domains[[ds]])
    domains[[de]] <- as.numeric(domains[[de]])
    ok_d <- is.finite(domains[[ds]]) & is.finite(domains[[de]]) &
        !is.na(domains[[dn]])
    if (any(!ok_d)) {
        message("Dropped ", sum(!ok_d), " domain row(s).")
        domains <- domains[ok_d, , drop = FALSE]
    }
    if (!nrow(domains)) {
        domains <- NULL
    } else {
        swap_d <- domains[[ds]] > domains[[de]]
        if (any(swap_d)) {
            tmp <- domains[[ds]][swap_d]
            domains[[ds]][swap_d] <- domains[[de]][swap_d]
            domains[[de]][swap_d] <- tmp
        }
        domains[[dn]] <- factor(
            as.character(domains[[dn]]),
            levels = unique(as.character(domains[[dn]]))
        )
        L <- max(L, max(domains[[de]], na.rm = TRUE))
    }
}

class_cols <- palette_colors("Qualitative.Bold", n = length(class_levels))
names(class_cols) <- class_levels

pop <- tolower(config$pop)
if (!pop %in% c("pie", "circle")) {
    stop("config$pop must be \"pie\" or \"circle\".", call. = FALSE)
}

max_n <- max(stems$total)
r_lo <- config$pop_r_mm[[1]]
r_hi <- config$pop_r_mm[[2]]
stems$r_mm <- r_lo + (r_hi - r_lo) * sqrt(stems$total / max_n)

agg$total <- stems$total[match(agg$pos, stems$pos)]
agg$r_mm <- stems$r_mm[match(agg$pos, stems$pos)]
agg <- agg[order(agg$pos, agg$class), ]
frac <- agg$n / agg$total
starts <- as.numeric(stats::ave(frac, agg$pos, FUN = function(z) {
    c(0, cumsum(z)[-length(z)])
}))
agg$from <- starts * 2 * pi
agg$to <- (starts + frac) * 2 * pi

maj <- agg[order(agg$pos, -agg$n, agg$class), ]
maj <- maj[!duplicated(maj$pos), c("pos", "class")]
stems$majority <- maj$class[match(stems$pos, maj$pos)]
stems$majority <- factor(as.character(stems$majority), levels = class_levels)

lab_df <- stems[!is.na(stems$label), , drop = FALSE]

x_limits <- c(0.5, L + 0.5)
y_top <- max_n * 1.22

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
out_w <- config$fig$width
out_h <- config$fig$height
if (identical(config$layout, "circular")) {
    cc <- config$circular
    if (identical(pop, "pie")) {
        message("layout = \"circular\" draws solid heads; pie heads are linear only.")
    }
    dir_sign <- if (isTRUE(cc$clockwise)) 1 else -1
    to_deg <- function(pos) cc$start_degree + dir_sign * 360 * (pos - 1) / L
    polar_xy <- function(r, deg) {
        data.frame(x = r * sin(deg * pi / 180), y = r * cos(deg * pi / 180))
    }
    r_in <- 1
    r_out <- 1 + cc$ring_width
    outward <- !identical(cc$stem_direction, "inward")
    stem_max <- if (outward) cc$stem_length else min(cc$stem_length, 0.7)
    to_r <- function(v) {
        if (outward) r_out + 0.04 + stem_max * v / max_n
        else r_in - 0.04 - stem_max * v / max_n
    }

    circle_path <- function(r) {
        cbind(polar_xy(r, seq(0, 360, length.out = 361)), r = r)
    }
    ring_vals <- if (is.null(cc$rings)) pretty(c(0, max_n), n = 3) else cc$rings
    ring_vals <- ring_vals[ring_vals > 0 & ring_vals <= max_n * 1.05]
    guide_df <- do.call(rbind, lapply(ring_vals, function(v) circle_path(to_r(v))))
    ring_lab <- cbind(polar_xy(to_r(ring_vals), cc$start_degree + 4), label = ring_vals)

    feat_poly <- NULL
    feat_key <- NULL
    if (!is.null(domains)) {
        ds <- config$domain$start
        de <- config$domain$end
        dn <- config$domain$name
        dg <- config$domain$group
        key <- if (!is.null(dg) && dg %in% names(domains)) domains[[dg]] else domains[[dn]]
        domains$key <- factor(as.character(key), levels = unique(as.character(key)))
        feat_poly <- do.call(rbind, lapply(seq_len(nrow(domains)), function(k) {
            arc <- seq(to_deg(domains[[ds]][[k]]), to_deg(domains[[de]][[k]] + 1),
                length.out = max(4, ceiling(80 * (domains[[de]][[k]] - domains[[ds]][[k]]) / L)))
            pts <- rbind(polar_xy(r_out, arc), polar_xy(r_in, rev(arc)))
            pts$id <- k
            pts$key <- domains$key[[k]]
            pts
        }))
        feat_key <- levels(domains$key)
    }

    stems$deg <- to_deg(stems$pos)
    stem_base <- polar_xy(if (outward) r_out else r_in, stems$deg)
    stem_tip <- polar_xy(to_r(stems$total), stems$deg)
    seg_df <- data.frame(x = stem_base$x, y = stem_base$y,
        xend = stem_tip$x, yend = stem_tip$y, majority = stems$majority)
    head_df <- cbind(stem_tip, majority = stems$majority)
    lab_c <- stems[!is.na(stems$label), , drop = FALSE]
    lab_head <- polar_xy(to_r(lab_c$total), lab_c$deg)
    lab_push <- polar_xy(if (outward) 0.16 else -0.16, lab_c$deg)
    lab_c <- cbind(lab_c, hx = lab_head$x, hy = lab_head$y,
        nx = lab_push$x, ny = lab_push$y)

    lim <- max(abs(c(guide_df$x, guide_df$y, head_df$x, head_df$y, r_out))) + 0.35
    p <- ggplot2::ggplot() +
        ggplot2::geom_path(data = guide_df, ggplot2::aes(x = x, y = y, group = r),
            colour = "grey78", linetype = "dashed", linewidth = 0.45)
    if (!is.null(feat_poly)) {
        dom_cols <- palette_colors("Qualitative.Bold", n = length(feat_key))
        names(dom_cols) <- feat_key
        p <- p +
            ggplot2::geom_polygon(data = feat_poly,
                ggplot2::aes(x = x, y = y, group = id, fill = key),
                colour = "grey15", linewidth = 0.25) +
            ggplot2::scale_fill_manual(values = dom_cols, name = NULL,
                guide = ggplot2::guide_legend(ncol = cc$legend_ncol, byrow = TRUE))
    }
    p <- p +
        ggplot2::geom_path(data = circle_path(r_in), ggplot2::aes(x = x, y = y),
            colour = "grey10", linewidth = 0.5) +
        ggplot2::geom_path(data = circle_path(r_out), ggplot2::aes(x = x, y = y),
            colour = "grey10", linewidth = 0.5) +
        ggplot2::geom_segment(data = seg_df,
            ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
            colour = "grey10", linewidth = 0.6, lineend = "butt") +
        ggplot2::geom_point(data = head_df,
            ggplot2::aes(x = x, y = y, colour = majority),
            size = cc$head_size, shape = 16) +
        ggplot2::scale_colour_manual(
            values = if (length(class_levels) > 1L) class_cols else
                setNames("grey10", class_levels),
            name = config$labels$fill,
            guide = if (length(class_levels) > 1L) "legend" else "none")
    p <- p +
        ggplot2::geom_text(data = ring_lab, ggplot2::aes(x = x, y = y, label = label),
            size = 2.6, colour = "grey55", hjust = 0, vjust = -0.3)
    if (nrow(lab_c)) {
        p <- p + ggrepel::geom_text_repel(
            data = lab_c, ggplot2::aes(x = hx, y = hy, label = label),
            nudge_x = lab_c$nx, nudge_y = lab_c$ny,
            size = cc$label_size, colour = "grey10",
            min.segment.length = Inf, box.padding = 0.15,
            point.size = cc$head_size, point.padding = 0.1,
            max.overlaps = Inf, seed = 1
        )
    }
    p <- p +
        ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
        ggplot2::labs(title = config$labels$title) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom",
            legend.text = ggplot2::element_text(size = 11),
            legend.key.size = ggplot2::unit(4, "mm"),
            plot.margin = ggplot2::margin(6, 6, 6, 6)
        )
    g <- p
    out_w <- cc$width
    out_h <- cc$height
} else {
    p_lolli <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = stems,
            ggplot2::aes(x = pos, xend = pos, y = 0, yend = total),
            colour = "grey40",
            linewidth = 0.4,
            lineend = "butt"
        )

    if (identical(pop, "pie")) {
        p_lolli <- p_lolli +
            geom_pie_mm(
                data = agg,
                ggplot2::aes(
                    x = pos,
                    y = total,
                    from = from,
                    to = to,
                    fill = class,
                    r_mm = r_mm
                ),
                colour = "white",
                linewidth = 0.15,
                inherit.aes = FALSE
            )
    } else {
        p_lolli <- p_lolli +
            ggplot2::geom_point(
                data = stems,
                ggplot2::aes(
                    x = pos,
                    y = total,
                    fill = majority,
                    size = total
                ),
                shape = 21,
                colour = "grey25",
                stroke = 0.35
            ) +
            ggplot2::scale_size_area(max_size = 5, guide = "none")
    }

    p_lolli <- p_lolli +
        ggplot2::scale_fill_manual(
            values = class_cols,
            name = config$labels$fill,
            drop = FALSE
        ) +
        ggplot2::scale_x_continuous(
            limits = x_limits,
            expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::coord_cartesian(ylim = c(0, y_top), clip = "off") +
        ggplot2::labs(
            title = config$labels$title,
            x = NULL,
            y = config$labels$y
        ) +
        ggplot2::guides(
            fill = ggplot2::guide_legend(
                override.aes = list(colour = "grey25", linewidth = 0.3, size = 4)
            )
        ) +
        ggprism::theme_prism() +
        ggplot2::theme(
            plot.background = ggplot2::element_rect(fill = NA, colour = NA),
            panel.background = ggplot2::element_rect(fill = NA, colour = NA),
            legend.background = ggplot2::element_blank(),
            legend.box.background = ggplot2::element_blank(),
            legend.key = ggplot2::element_blank(),
            legend.position = "right",
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(6, 8, 2, 8)
        )

    if (nrow(lab_df)) {
        p_lolli <- p_lolli +
            ggrepel::geom_text_repel(
                data = lab_df,
                ggplot2::aes(x = pos, y = total, label = label),
                size = 2.6,
                colour = "grey15",
                min.segment.length = 0,
                box.padding = 0.28,
                point.padding = 0.35,
                max.overlaps = Inf,
                segment.colour = "grey50",
                segment.size = 0.25,
                seed = 1
            )
    }

    backbone <- data.frame(start = 1, end = L)
    p_dom <- ggplot2::ggplot() +
        ggplot2::geom_rect(
            data = backbone,
            ggplot2::aes(xmin = start, xmax = end, ymin = 0.38, ymax = 0.62),
            fill = "grey78",
            colour = NA
        )

    if (!is.null(domains)) {
        ds <- config$domain$start
        de <- config$domain$end
        dn <- config$domain$name
        n_dom <- nlevels(domains[[dn]])
        dom_cols <- palette_colors("Brand.Algolia", n = n_dom)
        names(dom_cols) <- levels(domains[[dn]])
        domains$mid <- (domains[[ds]] + domains[[de]]) / 2
        p_dom <- p_dom +
            ggplot2::geom_rect(
                data = domains,
                ggplot2::aes(
                    xmin = .data[[ds]],
                    xmax = .data[[de]],
                    ymin = 0.12,
                    ymax = 0.88,
                    fill = .data[[dn]]
                ),
                colour = "grey20",
                linewidth = 0.25
            ) +
            ggplot2::geom_text(
                data = domains,
                ggplot2::aes(
                    x = mid,
                    y = 0.5,
                    label = .data[[dn]]
                ),
                size = 2.5,
                colour = "grey10",
                fontface = "bold"
            ) +
            ggplot2::scale_fill_manual(values = dom_cols, guide = "none")
    }

    p_dom <- p_dom +
        ggplot2::scale_x_continuous(
            limits = x_limits,
            expand = c(0, 0)
        ) +
        ggplot2::coord_cartesian(ylim = c(0, 1), clip = "off") +
        ggplot2::labs(x = config$labels$x, y = NULL) +
        ggprism::theme_prism() +
        ggplot2::theme(
            plot.background = ggplot2::element_rect(fill = NA, colour = NA),
            panel.background = ggplot2::element_rect(fill = NA, colour = NA),
            legend.position = "none",
            axis.line.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 8, 6, 8)
        )

    p <- p_lolli / p_dom + patchwork::plot_layout(heights = c(3.55, 0.72))
    g <- patchwork::patchworkGrob(p)
    bg <- which(g$layout$name == "background")
    if (length(bg)) {
        for (i in bg) {
            g$grobs[[i]]$gp$fill <- NA
            g$grobs[[i]]$gp$col <- NA
        }
    }
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(g, io$output, width = out_w, height = out_h)
