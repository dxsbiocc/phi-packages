#!/usr/bin/env Rscript

# Template-ID: scatter-umap-atlas
# Purpose: supplied embedding + aligned radial cluster hierarchy + annotation rings.
# Inputs: one row per cell (see CONFIG); clusters.tsv and hierarchy.tsv beside input.
# Output: PNG, PDF or SVG. Example: Rscript plot.R example.tsv output.pdf
# Class/Subtype: Rscript plot.R --class-subtype class_subtype/example.tsv output.png
# Dependencies: circlize, MASS, readr (and jsonlite via the common palette helper).
# Scientific assumptions: embedding, tree topology/heights, labels and scores are supplied.
# No Seurat integration, UMAP, clustering, marker selection, pseudoage or downloads.
# Geometry inspired by plot1cell and the user-supplied Nature atlas figure; not a wrapper.
# 中文：只做展示。树和年龄评分必须提前提供，不能从本图推断谱系或显著性。

local({
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (!length(file_arg)) stop("Run this file with Rscript.", call. = FALSE)
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

class_subtype <- identical(commandArgs(trailingOnly = TRUE)[1L], "--class-subtype")
io <- if (class_subtype) {
    args <- parse_io_args(n_input = 2L, input_names = c("--class-subtype", "cells"))
    list(input = args$input[[2L]], output = args$output)
} else {
    parse_io_args()
}

# CONFIG：适配数据时修改此处；全部颜色来自项目调色板。
config <- list(
    columns = list(id = "cell_id", x = "umap_1", y = "umap_2",
                   cluster = "cluster", celltype = "Celltype"),
    sidecars = list(clusters = "clusters.tsv", hierarchy = "hierarchy.tsv"),
    cluster_columns = list(id = "cluster", celltype = "Celltype", order = "order", count = "count"),
    tree_columns = list(node = "node", parent = "parent", height = "height"),
    categorical_rings = c("Region"),
    numeric_rings = c("Age", "CycleScore"),
    composition_rings = c("Group", "orig.ident"),
    insets = list(
        list(column = "Region", title = "Region", kind = "categorical"),
        list(column = "Marker", title = "Marker annotation", kind = "categorical"),
        list(column = "Pseudoage", title = "Pseudoage", kind = "numeric"),
        list(column = "CycleScore", title = "Cell cycle score", kind = "numeric")
    ),
    palette = "Qualitative.Prism",
    numeric_palettes = c(Age = "Quantitative.BluGrn", Pseudoage = "Quantitative.BluGrn",
                         CycleScore = "Gradient.YlOrRd"),
    # Age 与 Pseudoage 共用色阶，要求相同单位；设为 character() 可禁用。
    shared_age_columns = c("Age", "Pseudoage"),
    start_degree = 90,
    gap_degree = 2,
    end_gap_degree = 8,
    check_alignment = TRUE,
    max_alignment_error = 45, # degrees; this is a display-order check, not a test
    embedding_radius = 0.78,
    point_cex = 0.20,
    point_alpha = 0.50,
    # 与 plot1cell 一致：全部细胞共用密度面，相邻不同类型可连成同一轮廓。
    # scope="celltype" 仅用于分群轮廓；levels 可显式给出绝对密度阈值。
    contour = list(show = TRUE, scope = "global", grid = 200,
                   peak_fraction = c(0.15, 0.30), levels = NULL,
                   bandwidth_adjust = 1, line_palette = "Qualitative.Safe", line_index = 10),
    labels = list(show = TRUE, cex = 0.55),
    radii = list(tree_inner = 0.86, tree_outer = 0.95, ring_start = 0.956,
                 ring_height = 0.010, ring_gap = 0.004, count_height = 0.10),
    size = list(width = 10, height = 10, dpi = 300)
)

if (class_subtype) {
    config$columns$cluster <- config$cluster_columns$id <- "Subtype"
    config$columns$celltype <- config$cluster_columns$celltype <- "Class"
    config$categorical_rings <- config$numeric_rings <- config$composition_rings <- character()
    config$insets <- list()
    config$start_degree <- 263
    config$gap_degree <- 0
    config$end_gap_degree <- 14
    # The supplied leaf order follows the tree. The embedding retains its orientation.
    config$check_alignment <- FALSE
    config$embedding_radius <- 0.68
    config$point_cex <- 0.26
    config$point_alpha <- 0.82
    config$contour$scope <- "celltype"
    config$contour$peak_fraction <- 0.14
    config$labels$cex <- 0.84
    config$radii <- list(tree_inner = 0.73, tree_outer = 0.94, ring_start = 0.95,
                         ring_height = 0.074, ring_gap = 0.009, count_height = 0)
}

load_packages(c("circlize", "MASS", "readr"))

# DATA PREPARATION：严格核对输入，既不抽样也不补齐未知细胞。
nonempty <- function(x, name) {
    x <- as.character(x)
    if (any(is.na(x) | !nzchar(trimws(x)))) stop(paste(name, "must be non-empty."), call. = FALSE)
    x
}
finite_numeric <- function(x, name) {
    x <- suppressWarnings(as.numeric(x))
    if (any(!is.finite(x))) stop(paste(name, "must be finite numeric data."), call. = FALSE)
    x
}
read_sidecar <- function(name) {
    path <- file.path(dirname(io$input), name)
    if (!file.exists(path)) stop(paste("Missing required sidecar:", path), call. = FALSE)
    as.data.frame(read_table_auto(path), stringsAsFactors = FALSE)
}

cells <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
if (!nrow(cells)) stop("Cell table must not be empty.", call. = FALSE)
if (length(config$insets) > 4L) stop("At most four corner insets are supported.", call. = FALSE)
inset_names <- vapply(config$insets, function(x) x$column, character(1))
inset_kinds <- vapply(config$insets, function(x) x$kind, character(1))
if (any(!inset_kinds %in% c("categorical", "numeric"))) stop("Unknown inset kind.", call. = FALSE)
require_columns(cells, c(config$columns, config$composition_rings, inset_names))
cc <- config$columns
cells$.id <- nonempty(cells[[cc$id]], "cell_id")
if (anyDuplicated(cells$.id)) stop("cell_id must be unique.", call. = FALSE)
cells$.cluster <- nonempty(cells[[cc$cluster]], "cell cluster")
cells$.celltype <- nonempty(cells[[cc$celltype]], "cell Celltype")
cells$.x <- finite_numeric(cells[[cc$x]], "embedding x")
cells$.y <- finite_numeric(cells[[cc$y]], "embedding y")
cell_categorical <- unique(c(config$composition_rings, inset_names[inset_kinds == "categorical"]))
cell_numeric <- unique(inset_names[inset_kinds == "numeric"])
for (nm in cell_categorical) cells[[nm]] <- nonempty(cells[[nm]], paste("cell", nm))
for (nm in cell_numeric) cells[[nm]] <- finite_numeric(cells[[nm]], paste("cell", nm))

clusters <- read_sidecar(config$sidecars$clusters)
kc <- config$cluster_columns
require_columns(clusters, c(kc, config$categorical_rings, config$numeric_rings))
if (nrow(clusters) < 2L) stop("Need at least two clusters.", call. = FALSE)
clusters$.id <- nonempty(clusters[[kc$id]], "cluster id")
clusters$.celltype <- nonempty(clusters[[kc$celltype]], "cluster Celltype")
clusters$.order <- finite_numeric(clusters[[kc$order]], "cluster order")
clusters$.count <- finite_numeric(clusters[[kc$count]], "cluster count")
n_cluster <- nrow(clusters)
if (anyDuplicated(clusters$.id)) stop("Cluster IDs must be unique.", call. = FALSE)
if (!identical(sort(clusters$.order), as.numeric(seq_len(n_cluster)))) {
    stop("Cluster order must be a permutation of 1..N.", call. = FALSE)
}
clusters <- clusters[order(clusters$.order), , drop = FALSE]
if (!setequal(cells$.cluster, clusters$.id)) {
    stop("Cell clusters and clusters.tsv IDs must match exactly.", call. = FALSE)
}
observed <- as.numeric(table(factor(cells$.cluster, levels = clusters$.id)))
if (any(clusters$.count <= 0 | clusters$.count %% 1 != 0 | clusters$.count != observed)) {
    stop("Cluster count must equal the number of supplied cells in each cluster.", call. = FALSE)
}
if (any(cells$.celltype != clusters$.celltype[match(cells$.cluster, clusters$.id)])) {
    stop("Celltype disagrees between cell and cluster tables.", call. = FALSE)
}
if (anyDuplicated(rle(clusters$.celltype)$values)) {
    stop("Each Celltype must occupy one contiguous block in cluster order.", call. = FALSE)
}
for (nm in config$categorical_rings) clusters[[nm]] <- nonempty(clusters[[nm]], paste("cluster", nm))
for (nm in config$numeric_rings) clusters[[nm]] <- finite_numeric(clusters[[nm]], paste("cluster", nm))

tree <- read_sidecar(config$sidecars$hierarchy)
tc <- config$tree_columns
require_columns(tree, tc)
tree <- data.frame(node = nonempty(tree[[tc$node]], "tree node"),
                   parent = as.character(tree[[tc$parent]]),
                   height = finite_numeric(tree[[tc$height]], "tree height"))
tree$parent[is.na(tree$parent)] <- ""
if (anyDuplicated(tree$node)) stop("Tree node IDs must be unique.", call. = FALSE)
roots <- tree$node[!nzchar(tree$parent)]
if (length(roots) != 1L) stop("Hierarchy must have exactly one root.", call. = FALSE)
root <- roots[[1]]
if (any(!tree$parent[nzchar(tree$parent)] %in% tree$node)) stop("Unknown tree parent.", call. = FALSE)
parents <- setNames(tree$parent, tree$node)
heights <- setNames(tree$height, tree$node)
for (node in tree$node) {
    seen <- character()
    current <- node
    while (nzchar(current)) {
        if (current %in% seen) stop("Hierarchy contains a cycle.", call. = FALSE)
        seen <- c(seen, current)
        current <- parents[[current]]
    }
}
leaves <- setdiff(tree$node, tree$parent[nzchar(tree$parent)])
if (!setequal(leaves, clusters$.id)) stop("Tree leaves must match all cluster IDs exactly.", call. = FALSE)
nonroot <- nzchar(tree$parent)
if (any(heights < 0) || heights[[root]] <= 0 || any(heights[leaves] != 0) ||
    any(tree$height[nonroot] > heights[tree$parent[nonroot]])) {
    stop("Tree heights must be nonnegative, zero at leaves, and nonincreasing from parent to child.", call. = FALSE)
}
children <- split(tree$node[nonroot], tree$parent[nonroot])
leaf_indices <- local({
    cache <- new.env(parent = emptyenv())
    find <- function(node) {
        if (exists(node, envir = cache, inherits = FALSE)) return(get(node, envir = cache))
        result <- if (node %in% leaves) match(node, clusters$.id) else
            sort(unlist(lapply(children[[node]], find), use.names = FALSE))
        assign(node, result, envir = cache)
        result
    }
    setNames(lapply(tree$node, find), tree$node)
})
for (indices in leaf_indices) {
    if (max(indices) - min(indices) + 1L != length(indices)) {
        stop("Each subtree must occupy contiguous cluster order positions.", call. = FALSE)
    }
}

# 等比例放入中心圆；角度顺序检查不会修改原始坐标或旋转树。
cells$.px <- cells$.x - mean(range(cells$.x))
cells$.py <- cells$.y - mean(range(cells$.y))
extent <- max(sqrt(cells$.px^2 + cells$.py^2))
if (extent <= 0) stop("Embedding coordinates must not all be identical.", call. = FALSE)
cells$.px <- cells$.px * config$embedding_radius / extent
cells$.py <- cells$.py * config$embedding_radius / extent
gaps <- c(ifelse(head(clusters$.celltype, -1) != tail(clusters$.celltype, -1), config$gap_degree, 0),
          config$end_gap_degree)
if (sum(gaps) >= 360) stop("Configured gaps exhaust the circle.", call. = FALSE)
span <- (360 - sum(gaps)) / n_cluster
start <- config$start_degree - c(0, head(cumsum(span + gaps), -1))
angles <- data.frame(start = start * pi / 180, end = (start - span) * pi / 180,
                     mid = (start - span / 2) * pi / 180)
if (config$check_alignment) {
    cx <- vapply(clusters$.id, function(id) mean(cells$.px[cells$.cluster == id]), numeric(1))
    cy <- vapply(clusters$.id, function(id) mean(cells$.py[cells$.cluster == id]), numeric(1))
    delta <- atan2(cy, cx) - angles$mid
    error <- abs(atan2(sin(delta), cos(delta))) * 180 / pi
    if (any(error > config$max_alignment_error)) {
        stop("Cluster order does not align with the embedding polar angles; check order/start_degree or explicitly disable check_alignment.", call. = FALSE)
    }
}

# PALETTES：类别颜色按名称匹配；数值色阶在主图和四角共享。
celltypes <- unique(clusters$.celltype)
ct_colours <- setNames(palette_colors(config$palette, length(celltypes)), celltypes)
ink <- palette_colors("Brand.Emma")[[1]]
grey <- palette_colors("Qualitative.Safe")[[12]]
categorical_names <- unique(c(config$categorical_rings, cell_categorical))
cat_colours <- setNames(lapply(categorical_names, function(nm) {
    values <- unique(c(if (nm %in% names(clusters)) as.character(clusters[[nm]]),
                       if (nm %in% names(cells)) as.character(cells[[nm]])))
    setNames(palette_colors("Qualitative.Safe", length(values)), values)
}), categorical_names)
numeric_names <- unique(c(config$numeric_rings, cell_numeric))
numeric_limits <- setNames(lapply(numeric_names, function(nm) {
    range(c(if (nm %in% config$numeric_rings) clusters[[nm]],
            if (nm %in% cell_numeric) cells[[nm]]))
}), numeric_names)
shared <- intersect(config$shared_age_columns, numeric_names)
if (length(shared) > 1L) {
    limits <- range(unlist(numeric_limits[shared], use.names = FALSE))
    for (nm in shared) numeric_limits[[nm]] <- limits
}
numeric_colours <- setNames(lapply(numeric_names, function(nm) {
    limits <- numeric_limits[[nm]]
    # 常量只扩展色阶，不改变数据值。
    if (diff(limits) == 0) limits <- limits + c(-1, 1) * max(abs(limits[[1]]) * 0.01, 0.01)
    if (!nm %in% names(config$numeric_palettes)) {
        stop(paste("Set numeric_palettes for", nm), call. = FALSE)
    }
    palette <- config$numeric_palettes[[nm]]
    cols <- palette_colors(palette)
    circlize::colorRamp2(seq(limits[[1]], limits[[2]], length.out = length(cols)), cols)
}), numeric_names)
n_ring <- if (class_subtype) 2L else
    1L + length(config$categorical_rings) + length(config$numeric_rings) + length(config$composition_rings)
ring_radii <- config$radii$ring_start + (seq_len(n_ring) - 1) *
    (config$radii$ring_height + config$radii$ring_gap)
count_base <- tail(ring_radii, 1) + config$radii$ring_height + 0.012
if (count_base + config$radii$count_height > 1.18) stop("Too many rings; reduce CONFIG ring heights/gaps.", call. = FALSE)

# PLOT HELPERS：全部环、树和小图都是可编辑的矢量元素。
polar <- function(theta, radius) cbind(cos(theta) * radius, sin(theta) * radius)
annulus <- function(a, b, r0, r1, colour) {
    theta <- seq(a, b, length.out = max(3L, ceiling(abs(a - b) * 120)))
    xy <- rbind(polar(theta, r0), polar(rev(theta), r1))
    graphics::polygon(xy, col = colour, border = NA)
}
draw_tree <- function() {
    radius <- config$radii$tree_inner +
        (1 - heights / heights[[root]]) * (config$radii$tree_outer - config$radii$tree_inner)
    theta <- vapply(leaf_indices, function(ids) mean(range(angles$mid[ids])), numeric(1))
    branch_colours <- if (class_subtype) vapply(leaf_indices, function(ids) {
        classes <- unique(clusters$.celltype[ids])
        if (length(classes) == 1L) ct_colours[[classes]] else ink
    }, character(1)) else setNames(rep(grey, length(theta)), names(theta))
    branch_width <- if (class_subtype) 1.1 else 0.45
    for (node in names(children)) {
        kids <- children[[node]]
        child_angles <- theta[kids]
        arc <- seq(min(child_angles), max(child_angles), length.out = 60)
        graphics::lines(polar(arc, radius[[node]]), col = branch_colours[[node]], lwd = branch_width)
        for (child in kids) {
            xy <- polar(rep(theta[[child]], 2), c(radius[[node]], radius[[child]]))
            graphics::lines(xy, col = branch_colours[[child]], lwd = branch_width)
        }
    }
}
draw_rings <- function() {
    ring <- 0L
    add <- function(colours) {
        ring <<- ring + 1L
        r <- ring_radii[[ring]]
        for (j in seq_len(n_cluster)) annulus(angles$start[[j]], angles$end[[j]],
            r, r + config$radii$ring_height, colours[[j]])
    }
    add(unname(ct_colours[clusters$.celltype]))
    for (nm in config$categorical_rings) add(unname(cat_colours[[nm]][clusters[[nm]]]))
    for (nm in config$numeric_rings) add(numeric_colours[[nm]](clusters[[nm]]))
    for (nm in config$composition_rings) {
        ring <- ring + 1L
        r <- ring_radii[[ring]]
        levels <- names(cat_colours[[nm]])
        for (j in seq_len(n_cluster)) {
            counts <- as.numeric(table(factor(cells[[nm]][cells$.cluster == clusters$.id[[j]]], levels = levels)))
            ends <- cumsum(counts / sum(counts))
            starts <- c(0, head(ends, -1))
            for (k in which(counts > 0)) annulus(
                angles$start[[j]] - starts[[k]] * (angles$start[[j]] - angles$end[[j]]),
                angles$start[[j]] - ends[[k]] * (angles$start[[j]] - angles$end[[j]]),
                r, r + config$radii$ring_height, cat_colours[[nm]][[k]])
        }
    }
    for (j in seq_len(n_cluster)) {
        width <- (angles$start[[j]] - angles$end[[j]]) * 0.28
        annulus(angles$mid[[j]] + width, angles$mid[[j]] - width, count_base,
            count_base + config$radii$count_height * clusters$.count[[j]] / max(clusters$.count), grey)
    }
    for (type in celltypes) {
        ids <- which(clusters$.celltype == type)
        theta <- mean(range(angles$mid[ids]))
        xy <- polar(theta, count_base + config$radii$count_height + 0.045)
        rotation <- (theta * 180 / pi) %% 180 - 90
        graphics::text(xy, labels = type, srt = rotation, cex = 0.53, col = ink)
    }
}
draw_class_subtype_rings <- function() {
    height <- config$radii$ring_height
    class_radius <- ring_radii[[1L]]
    subtype_radius <- ring_radii[[2L]]
    for (type in celltypes) {
        ids <- which(clusters$.celltype == type)
        a <- angles$start[min(ids)]
        b <- angles$end[max(ids)]
        colour <- ct_colours[[type]]
        annulus(a, b, class_radius, class_radius + height, colour)
        theta <- (a + b) / 2
        label_colour <- if (.hex_luminance(colour) < 0.48) "white" else ink
        label_cex <- min(0.68, 0.90 * (a - b) * class_radius / graphics::strwidth(type))
        graphics::text(polar(theta, class_radius + height / 2), labels = type,
                       srt = (theta * 180 / pi) %% 180 - 90, cex = label_cex, col = label_colour)
    }
    subtype_colours <- palette_colors("Qualitative.Safe", n_cluster)
    # Shrink the labels as the available arc per subtype gets shorter.
    label_radius <- subtype_radius + height + 0.022
    subtype_cex <- min(0.62, 15 / n_cluster,
                       (1.46 - label_radius) / max(graphics::strwidth(clusters$.id)))
    for (j in seq_len(n_cluster)) {
        annulus(angles$start[j], angles$end[j], subtype_radius, subtype_radius + height,
                subtype_colours[j])
        theta <- angles$mid[j]
        rotation <- theta * 180 / pi
        # Anchor labels at the rim and read upright on both sides of the circle.
        flip <- cos(theta) < 0
        graphics::text(polar(theta, label_radius), labels = clusters$.id[j],
                       srt = rotation + if (flip) 180 else 0, adj = c(as.numeric(flip), 0.5),
                       cex = subtype_cex, col = ink)
    }
    gap_mid <- tail(angles$end, 1L) - config$end_gap_degree * pi / 360
    graphics::text(polar(gap_mid, ring_radii + height / 2), labels = c("Class", "Subtype"),
                   cex = 0.52, col = ink)
}

draw_class_contours <- function(paths, colour) {
    # Fill each density level independently. The even-odd rule preserves holes.
    levels <- split(paths, vapply(paths, function(path) path$level, numeric(1)))
    for (level in levels) {
        closed <- vapply(level, function(path) {
            isTRUE(all.equal(c(path$x[1L], path$y[1L]), c(tail(path$x, 1L), tail(path$y, 1L))))
        }, logical(1))
        if (any(closed)) {
            x <- unlist(lapply(level[closed], function(path) c(path$x, NA_real_)))
            y <- unlist(lapply(level[closed], function(path) c(path$y, NA_real_)))
            graphics::polypath(head(x, -1L), head(y, -1L), rule = "evenodd",
                               col = grDevices::adjustcolor(colour, alpha.f = 0.50),
                               border = ink, lwd = 0.85)
        }
        for (path in level[!closed]) graphics::lines(path$x, path$y, col = ink, lwd = 0.85)
    }
}

build_embedding_contours <- function() {
    settings <- config$contour
    if (!settings$scope %in% c("global", "celltype")) {
        stop("contour scope must be global or celltype.", call. = FALSE)
    }
    if (length(settings$grid) != 1L || !is.finite(settings$grid) ||
        settings$grid < 20 || settings$grid %% 1 != 0) {
        stop("contour grid must be an integer >= 20.", call. = FALSE)
    }
    if (length(settings$bandwidth_adjust) != 1L || !is.finite(settings$bandwidth_adjust) ||
        settings$bandwidth_adjust <= 0) {
        stop("contour bandwidth_adjust must be positive and finite.", call. = FALSE)
    }
    requested <- if (is.null(settings$levels)) settings$peak_fraction else settings$levels
    if (!length(requested) || any(!is.finite(requested) | requested <= 0) ||
        (is.null(settings$levels) && any(requested >= 1))) {
        stop("Contour levels must be positive; peak_fraction must be between 0 and 1.", call. = FALSE)
    }
    groups <- if (settings$scope == "global") list("All cells" = cells) else
        split(cells, cells$.celltype, drop = TRUE)
    setNames(lapply(names(groups), function(name) {
        d <- groups[[name]]
        if (nrow(d) < 3L || length(unique(d$.px)) < 2L || length(unique(d$.py)) < 2L) {
            message("Skipped display contour for ", name, ": too few varying coordinates; all cells retained.")
            return(list(paths = list(), n_cells = nrow(d)))
        }
        bandwidth <- c(MASS::bandwidth.nrd(d$.px), MASS::bandwidth.nrd(d$.py)) * settings$bandwidth_adjust
        if (any(!is.finite(bandwidth) | bandwidth <= 0)) {
            message("Skipped display contour for ", name, ": non-positive bandwidth; all cells retained.")
            return(list(paths = list(), n_cells = nrow(d)))
        }
        # MASS::kde2d divides the supplied bandwidth by four internally.
        # Padding the grid by three actual bandwidths avoids cropped contour paths.
        limit <- config$embedding_radius + 3 * max(bandwidth) / 4
        kde <- MASS::kde2d(d$.px, d$.py, n = settings$grid, h = bandwidth,
                          lims = c(-limit, limit, -limit, limit))
        levels <- if (is.null(settings$levels)) max(kde$z) * requested else requested
        if (any(levels >= max(kde$z))) {
            stop(paste("Requested contour level is not below peak density for", name), call. = FALSE)
        }
        list(paths = grDevices::contourLines(kde$x, kde$y, kde$z, levels = sort(unique(levels))),
             n_cells = nrow(d), levels = levels, bandwidth = bandwidth)
    }), names(groups))
}

draw_embedding <- function() {
    contours <- if (config$contour$show) build_embedding_contours() else list()
    if (config$contour$show && config$contour$scope == "celltype") {
        for (type in names(contours)) {
            if (class_subtype) {
                draw_class_contours(contours[[type]]$paths, ct_colours[[type]])
                next
            }
            for (line in contours[[type]]$paths) {
                graphics::polygon(line$x, line$y,
                    col = grDevices::adjustcolor(ct_colours[[type]], alpha.f = 0.10),
                    border = grDevices::adjustcolor(ct_colours[[type]], alpha.f = 0.85), lwd = 0.65)
            }
        }
    }
    graphics::points(cells$.px, cells$.py, pch = 16, cex = config$point_cex,
        col = grDevices::adjustcolor(unname(ct_colours[cells$.celltype]), alpha.f = config$point_alpha))
    if (config$contour$show && config$contour$scope == "global") {
        colour <- palette_colors(config$contour$line_palette)[[config$contour$line_index]]
        # Shared outlines follow the pooled density, not Celltype boundaries.
        # Draw paths without polygon fill so holes in the density remain holes.
        for (line in contours[[1]]$paths) {
            line_width <- if (line$level == min(contours[[1]]$levels)) 1.1 else 0.8
            graphics::lines(line$x, line$y, col = grDevices::adjustcolor(colour, alpha.f = 0.85),
                            lwd = line_width)
        }
    }
    if (config$labels$show) {
        for (type in celltypes) {
            d <- cells[cells$.celltype == type, , drop = FALSE]
            x <- median(d$.px)
            y <- median(d$.py)
            label <- paste(strwrap(type, width = 14), collapse = "\n")
            for (offset in list(c(-0.004, 0), c(0.004, 0), c(0, -0.004), c(0, 0.004))) {
                graphics::text(x + offset[[1]], y + offset[[2]], label,
                               cex = config$labels$cex, col = "white")
            }
            graphics::text(x, y, label, cex = config$labels$cex, col = ink,
                           font = if (class_subtype) 2 else 1)
        }
    }
    if (!class_subtype) {
        graphics::text(0, config$radii$tree_inner - 0.035, paste0(n_cluster, " clusters  |  ",
            format(nrow(cells), big.mark = ","), " cells"), cex = 0.60, col = ink)
    }
}
draw_inset <- function(spec, position) {
    graphics::par(fig = position, mar = c(0, 0, 0, 0))
    # par(bg=...) resets new on PDF devices; set overlay state last.
    graphics::par(new = TRUE)
    graphics::plot.new()
    graphics::plot.window(xlim = c(-0.78, 0.78), ylim = c(-1.04, 0.88), asp = 1)
    nm <- spec$column
    colours <- if (spec$kind == "numeric") numeric_colours[[nm]](cells[[nm]]) else
        unname(cat_colours[[nm]][cells[[nm]]])
    graphics::points(cells$.px * 0.62 / config$embedding_radius,
                     cells$.py * 0.62 / config$embedding_radius, pch = 16, cex = 0.10,
                     col = grDevices::adjustcolor(colours, alpha.f = 0.6))
    graphics::text(0, 0.80, spec$title, cex = 0.64, col = ink, font = 2)
    if (spec$kind == "categorical") {
        map <- cat_colours[[nm]]
        graphics::legend("bottomleft", legend = names(map), col = unname(map),
            pch = 16, pt.cex = 0.6, cex = 0.46, ncol = 2, bty = "n", text.col = ink,
            x.intersp = 0.4, y.intersp = 0.85)
    } else {
        limits <- numeric_limits[[nm]]
        x <- seq(-0.60, 0.60, length.out = 41)
        values <- seq(limits[[1]], limits[[2]], length.out = 40)
        graphics::rect(head(x, -1), -0.87, tail(x, -1), -0.80,
                        col = numeric_colours[[nm]](values), border = NA)
        graphics::text(c(-0.60, 0, 0.60), -0.96,
            format(signif(c(limits[[1]], mean(limits), limits[[2]]), 3), trim = TRUE),
            cex = 0.46, col = ink)
    }
}
draw_atlas <- function() {
    graphics::par(mar = c(0.1, 0.1, 0.1, 0.1), bg = "transparent", xaxs = "i", yaxs = "i")
    graphics::plot.new()
    graphics::plot.window(xlim = c(-1.50, 1.50), ylim = c(-1.50, 1.50), asp = 1)
    draw_tree()
    if (class_subtype) {
        draw_class_subtype_rings()
        draw_embedding()
        return(invisible(NULL))
    }
    draw_rings()
    draw_embedding()
    key <- c("Tree", "Class", config$categorical_rings, config$numeric_rings,
             config$composition_rings, "Cluster size")
    graphics::text(0, -1.25, paste("Inner -> outer:", paste(key, collapse = " | ")),
                   cex = 0.44, col = ink)
    graphics::text(0, -1.30, paste0("Bar scale: 0-", max(clusters$.count), " cells per cluster"),
                   cex = 0.46, col = ink)
    y <- -1.36
    for (nm in config$composition_rings) {
        map <- cat_colours[[nm]]
        graphics::legend(x = 0, y = y, legend = c(paste0(nm, ":"), names(map)),
            col = c(NA, unname(map)), pch = c(NA, rep(15, length(map))),
            horiz = TRUE, xjust = 0.5, cex = 0.40, pt.cex = 0.6,
            bty = "n", text.col = ink, x.intersp = 0.4)
        y <- y - 0.07
    }
    positions <- list(c(0.015, 0.205, 0.785, 0.995), c(0.795, 0.985, 0.785, 0.995),
                       c(0.015, 0.205, 0.005, 0.215), c(0.795, 0.985, 0.005, 0.215))
    for (i in seq_along(config$insets)) draw_inset(config$insets[[i]], positions[[i]])
}

# SAVE：保留矢量 PDF/SVG，并输出 300 dpi 透明 PNG 预览。
fmt <- tolower(tools::file_ext(io$output))
if (!fmt %in% c("png", "pdf", "svg")) stop("Output must be png, pdf or svg.", call. = FALSE)
dir.create(dirname(io$output), recursive = TRUE, showWarnings = FALSE)
if (fmt == "png") {
    grDevices::png(io$output, width = config$size$width, height = config$size$height,
                   units = "in", res = config$size$dpi, bg = "transparent")
} else if (fmt == "pdf") {
    grDevices::pdf(io$output, width = config$size$width, height = config$size$height, bg = "transparent")
} else {
    grDevices::svg(io$output, width = config$size$width, height = config$size$height, bg = "transparent")
}
tryCatch(draw_atlas(), finally = grDevices::dev.off())
message("Saved: ", normalizePath(io$output, mustWork = TRUE))
