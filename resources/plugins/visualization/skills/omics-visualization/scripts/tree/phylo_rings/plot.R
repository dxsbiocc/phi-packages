#!/usr/bin/env Rscript

# Template-ID: tree-phylo-rings
#
# Purpose:
#   iTOL-style multi-ring phylogeny: a supplied Newick tree in an open
#   fan layout with clade background colour by a supplied group (e.g.
#   family), tip labels, one or more heatmap rings of supplied tip
#   values, and one or more stacked-bar rings of supplied per-tip
#   category values with dotted radial gridlines and tick labels in
#   the opening gap.
#
# Inputs:
#   tree.nwk - Newick tree (tip labels must match the tables)
#   tips.tsv - one row per tip: label, group, plus heatmap columns
#   bars.tsv - long table: label, track, category, value
#              (header-only = no bar rings)
#
# Output:
#   A PDF, PNG, or SVG circular phylogeny with annotation rings.
#
# Dependencies:
#   ggplot2, readr, ape, ggtree, ggnewscale, ComplexHeatmap, circlize, grid
#
# Example:
#   Rscript plot.R tree.nwk tips.tsv bars.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, heatmap rings and their transform /
#   palette, bar tracks and their palettes, open angle, ring widths,
#   label band). Extra rings are extra CONFIG entries; extra tips are
#   extra rows. Do not build, root, or prune the tree, or count genomes
#   / BGCs in this script.
#
# Scientific assumptions:
#   Tree topology, branch lengths, tip groups and all ring values are
#   supplied. Heatmap transform (e.g. log2(x + 1)) is a display
#   transform. Clade shading uses the MRCA of each group's tips; a
#   non-monophyletic group is shaded tip by tip instead.

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

io <- parse_io_args(n_input = 3L, input_names = c("tree", "tips", "bars"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    tip_columns = list(label = "label", group = "family"),
    group_title = "Family",
    group_palette = "Qualitative.Pastel",
    group_alpha = 0.75,
    # heatmap rings, inner -> outer
    heat = list(
        list(column = "genome_count", title = "Genome count (log2)",
            transform = "log2p1", palette = "Gradient.Greys"),
        list(column = "bgc_count", title = "BGC count (log2)",
            transform = "log2p1", palette = "Gradient.RdPu")
    ),
    bar_columns = list(label = "label", track = "track",
        category = "category", value = "value"),
    # stacked-bar rings, inner -> outer; NULL palette = Qualitative.Bold
    bars = list(
        list(track = "Average BGC count per genome",
            palette = "Qualitative.Bold"),
        list(track = "Average RiPP count per genome",
            palette = "Qualitative.Prism")
    ),
    underscore_to_space = TRUE,  # Newick writes spaces as "_"
    open_angle = 90,             # leave first quadrant empty for the legend
    rotate = 90,                 # shift the open sector into Q1 (top-right)
    label_size = 2.4,
    label_band = 0.55,           # tip-label band, fraction of tree depth
    heat_width = 0.07,           # each heatmap ring, fraction of tree depth
    bar_width = 0.32,            # each bar ring, fraction of tree depth
    ring_gap = 0.02,
    # Legend sits in the open first quadrant as a 2-column rectangle:
    # left = Family + heatmap colourbars; right = stacked-bar legends.
    legend = list(
        x = 0.72, y = 0.72,      # npc centre of the open first quadrant
        width = 0.36, height = 0.36
    ),
    size = list(
        width = 10,
        height = 10
    )
)

load_packages(c("ggplot2", "readr", "ape", "ggtree", "ggnewscale",
    "ComplexHeatmap", "circlize", "grid"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
tr <- ape::read.tree(io$input[[1]])
if (is.null(tr$edge.length)) stop("Tree needs branch lengths.", call. = FALSE)
tips <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
tc <- config$tip_columns
heat_cols <- vapply(config$heat, `[[`, "", "column")
require_columns(tips, c(tc, as.list(heat_cols)))
tips$label <- as.character(tips[[tc$label]])
if (isTRUE(config$underscore_to_space)) {
    tr$tip.label <- gsub("_", " ", tr$tip.label, fixed = TRUE)
    tips$label <- gsub("_", " ", tips$label, fixed = TRUE)
}
tips$group <- as.character(tips[[tc$group]])
miss <- setdiff(tr$tip.label, tips$label)
if (length(miss)) {
    stop(sprintf("Tips missing from tips table: %s",
        paste(utils::head(miss, 5), collapse = ", ")), call. = FALSE)
}
tips <- tips[match(tr$tip.label, tips$label), , drop = FALSE]

bars <- as.data.frame(read_table_auto(io$input[[3]]), stringsAsFactors = FALSE)
bc <- config$bar_columns
if (nrow(bars)) {
    require_columns(bars, bc)
    bars <- data.frame(label = as.character(bars[[bc$label]]),
        track = as.character(bars[[bc$track]]),
        category = as.character(bars[[bc$category]]),
        value = as.numeric(bars[[bc$value]]), stringsAsFactors = FALSE)
    if (any(!is.finite(bars$value) | bars$value < 0)) {
        stop("Bar values must be finite and non-negative.", call. = FALSE)
    }
    if (isTRUE(config$underscore_to_space)) {
        bars$label <- gsub("_", " ", bars$label, fixed = TRUE)
    }
    bars <- bars[bars$label %in% tr$tip.label, , drop = FALSE]
}

base <- ggtree::ggtree(tr, layout = "fan", open.angle = config$open_angle,
    linewidth = 0.3)
td <- base$data
tip_y <- setNames(td$y[td$isTip], td$label[td$isTip])
depth <- max(td$x)
n_tip <- length(tip_y)

# radii of every ring
r_label_end <- depth * (1 + config$label_band)
gap <- depth * config$ring_gap
hw <- depth * config$heat_width
heat_x <- r_label_end + gap + (seq_along(config$heat) - 0.5) * (hw + gap)
bar_base <- r_label_end + length(config$heat) * (hw + gap) + gap * 2
bw <- depth * config$bar_width
bar_x0 <- bar_base + (seq_along(config$bars) - 1) * (bw + gap * 3)

grp_lv <- unique(tips$group)
grp_col <- setNames(palette_colors(config$group_palette, n = length(grp_lv)),
    grp_lv)

# clade shading: MRCA node when monophyletic, otherwise per-tip rects
hl <- list()
tip_rects <- list()
for (g in grp_lv) {
    labs_g <- tips$label[tips$group == g]
    node <- if (length(labs_g) > 1) ape::getMRCA(tr, labs_g) else
        which(tr$tip.label == labs_g)
    desc <- if (node <= ape::Ntip(tr)) tr$tip.label[node] else
        ape::extract.clade(tr, node)$tip.label
    if (length(labs_g) > 1 && setequal(desc, labs_g)) {
        hl[[g]] <- data.frame(node = node, group = g)
    } else {
        tip_rects[[g]] <- data.frame(group = g, y = tip_y[labs_g])
    }
}
hl <- do.call(rbind, hl)
tip_rects <- do.call(rbind, tip_rects)

heat_df <- lapply(seq_along(config$heat), function(i) {
    h <- config$heat[[i]]
    v <- as.numeric(tips[[h$column]])
    if (identical(h$transform, "log2p1")) v <- log2(v + 1)
    data.frame(y = tip_y[tips$label], x = heat_x[[i]], value = v)
})

bar_df <- lapply(seq_along(config$bars), function(i) {
    b <- config$bars[[i]]
    d <- bars[bars$track == b$track, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    cats <- unique(d$category)
    d$category <- factor(d$category, levels = cats)
    d <- d[order(d$label, d$category), , drop = FALSE]
    d$hi <- stats::ave(d$value, d$label, FUN = cumsum)
    d$lo <- d$hi - d$value
    tot <- max(tapply(d$value, d$label, sum))
    tick <- pretty(c(0, tot), n = 3)
    tick <- tick[tick <= tot * 1.0001]
    sc <- bw / tot
    d$xmin <- bar_x0[[i]] + d$lo * sc
    d$xmax <- bar_x0[[i]] + d$hi * sc
    d$y <- tip_y[d$label]
    list(data = d, ticks = data.frame(r = bar_x0[[i]] + tick * sc, lab = tick),
        cats = cats)
})

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- base
if (!is.null(hl)) {
    p <- p + ggtree::geom_hilight(data = hl, ggplot2::aes(node = node,
        fill = group), extendto = r_label_end, alpha = config$group_alpha,
        to.bottom = TRUE)
}
if (!is.null(tip_rects)) {
    p <- p + geom_tile(data = tip_rects, aes(x = depth * (1 + config$label_band / 2),
        y = y, fill = group), width = depth * config$label_band, height = 1,
        alpha = config$group_alpha, inherit.aes = FALSE)
}
p <- p +
    scale_fill_manual(values = grp_col, guide = "none") +
    ggtree::geom_tiplab(size = config$label_size, offset = depth * 0.02)

heat_funs <- list()
for (i in seq_along(heat_df)) {
    h <- config$heat[[i]]
    stops <- palette_colors(h$palette)
    lims <- range(heat_df[[i]]$value, na.rm = TRUE)
    heat_funs[[i]] <- list(
        title = h$title,
        col_fun = circlize::colorRamp2(
            seq(lims[[1]], lims[[2]], length.out = length(stops)), stops),
        at = pretty(lims, n = 4)
    )
    p <- p + ggnewscale::new_scale_fill() +
        geom_tile(data = heat_df[[i]], aes(x = x, y = y, fill = value),
            width = hw, height = 0.92, colour = "white", linewidth = 0.15,
            inherit.aes = FALSE) +
        scale_fill_gradientn(colours = stops, guide = "none")
}

bar_legs <- list()
for (i in seq_along(bar_df)) {
    bd <- bar_df[[i]]
    if (is.null(bd)) next
    pal_id <- config$bars[[i]]$palette
    if (is.null(pal_id)) pal_id <- "Qualitative.Bold"
    cc <- setNames(palette_colors(pal_id, n = length(bd$cats)), bd$cats)
    bar_legs[[i]] <- list(title = config$bars[[i]]$track, colors = cc)
    p <- p +
        geom_vline(xintercept = bd$ticks$r, linetype = "dotted",
            colour = "grey70", linewidth = 0.25) +
        ggnewscale::new_scale_fill() +
        geom_rect(data = bd$data, aes(xmin = xmin, xmax = xmax,
            ymin = y - 0.42, ymax = y + 0.42, fill = category),
            inherit.aes = FALSE) +
        scale_fill_manual(values = cc, guide = "none") +
        annotate("text", x = bd$ticks$r, y = n_tip + 1.2,
            label = format(bd$ticks$lab), size = 2.6, angle = 0)
}

# Keep fan open.angle; only rotate the polar origin so the empty
# sector lands in the first quadrant (top-right). Tip-label angles
# track the same offset so text stays radial.
p$data$angle <- p$data$angle + config$rotate
p <- p +
    ggplot2::coord_polar(theta = "y",
        start = (-90 + config$rotate) / 180 * pi, direction = -1) +
    theme(
        legend.position = "none",
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(4, 4, 4, 4)
    )

build_legend <- function() {
    left <- list(ComplexHeatmap::Legend(
        labels = names(grp_col),
        legend_gp = grid::gpar(fill = unname(grp_col)),
        title = config$group_title,
        title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 8),
        grid_height = grid::unit(3.2, "mm"),
        grid_width = grid::unit(3.2, "mm"),
        row_gap = grid::unit(0.8, "mm")
    ))
    for (hf in heat_funs) {
        left[[length(left) + 1L]] <- ComplexHeatmap::Legend(
            col_fun = hf$col_fun, at = hf$at, title = hf$title,
            title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = 7.5),
            direction = "horizontal",
            legend_width = grid::unit(28, "mm"),
            grid_height = grid::unit(3, "mm")
        )
    }
    right <- list()
    for (bl in bar_legs) {
        if (is.null(bl)) next
        right[[length(right) + 1L]] <- ComplexHeatmap::Legend(
            labels = names(bl$colors),
            legend_gp = grid::gpar(fill = unname(bl$colors)),
            title = bl$title,
            title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = 8),
            grid_height = grid::unit(3.2, "mm"),
            grid_width = grid::unit(3.2, "mm"),
            row_gap = grid::unit(0.8, "mm")
        )
    }
    left_pack <- do.call(ComplexHeatmap::packLegend,
        c(left, list(direction = "vertical", gap = grid::unit(3, "mm"))))
    if (length(right)) {
        right_pack <- do.call(ComplexHeatmap::packLegend,
            c(right, list(direction = "vertical", gap = grid::unit(3, "mm"))))
        ComplexHeatmap::packLegend(left_pack, right_pack,
            direction = "horizontal", gap = grid::unit(5, "mm"))
    } else {
        left_pack
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
        res = 300, bg = "white")
} else if (identical(fmt, "pdf")) {
    grDevices::pdf(out_path, width = wd, height = ht)
} else {
    grDevices::svg(out_path, width = wd, height = ht, bg = "white")
}
tryCatch({
    print(p)
    lg <- config$legend
    grid::pushViewport(grid::viewport(x = lg$x, y = lg$y,
        width = lg$width, height = lg$height, just = c(0.5, 0.5)))
    ComplexHeatmap::draw(build_legend(), x = grid::unit(0.5, "npc"),
        y = grid::unit(0.5, "npc"), just = c("center", "center"))
    grid::popViewport()
}, finally = grDevices::dev.off())
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
