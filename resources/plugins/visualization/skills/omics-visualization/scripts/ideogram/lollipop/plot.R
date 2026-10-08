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
#     - domains.tsv: start, end, domain
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
        name = "domain"
    ),
    protein_length = 393,
    pop = "pie",
    pop_r_mm = c(1.0, 3.0),
    label_min = 8,
    fig = list(
        width = 8.8,
        height = 4.2
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
stems$label[need_lab] <- paste0("p.", stems$pos[need_lab])
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
    require_columns(domains, config$domain)
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

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(g, io$output, width = config$fig$width, height = config$fig$height)
