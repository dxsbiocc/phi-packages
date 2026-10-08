#!/usr/bin/env Rscript

# Template-ID: heatmap-cluster-block
#
# Purpose:
#   Draw a clustered rectangular heatmap with dendrograms cut into
#   k row and column blocks (pheatmap-style cutree gaps) and an
#   optional row-group annotation strip.
#
# Inputs:
#   A table with one row-id column and numeric score columns
#   (typically a supplied Pearson matrix: cell type × factor).
#   Default example:
#     - cell_type: row label
#     - remaining columns: one numeric score per column id
#   Optional table beside the input (same directory):
#     - rowInfo.tsv: row id, Group
#
# Output:
#   A PDF, PNG, or SVG clustered heatmap with block gaps.
#
# Dependencies:
#   ComplexHeatmap, circlize, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (id column, cutree k, labels,
#   fill limits, group colours). Place rowInfo.tsv next to the input;
#   the CLI stays one input plus one output. Edit DATA PREPARATION
#   to change distance or linkage. Edit PLOT only when the split
#   geometry must change. cutree k is CONFIG, not a new id.
#
# Scientific assumptions:
#   Each cell is a supplied number. This script does not compute
#   correlation, NMF, or z-score. Dendrograms and cutree blocks are
#   display clustering of the supplied matrix, not an upstream test.

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
        id = "cell_type"
    ),
    cutree_rows = 3,
    cutree_cols = 3,
    clustering_method = "complete",
    fill_limits = c(-0.5, 0.5),
    labels = list(
        title = "NMF factors",
        y = "Mapped cell type",
        fill = "Pearson cor.",
        group = "Group"
    ),
    # Qualitative.Prism (orange, purple, yellow, blue, green)
    group_colors = c(
        Epithelium = "#e17c05",
        Endothelium = "#5f4690",
        Mesenchyme = "#edad08",
        Lymphoid = "#1d6996",
        Myeloid = "#73af48"
    )
)

load_packages(c("ComplexHeatmap", "circlize", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

id_col <- config$columns$id
value_cols <- setdiff(names(df), id_col)
if (length(value_cols) < 2L) {
    stop("Heatmap input needs at least two numeric columns besides the id.",
        call. = FALSE
    )
}

ok <- !is.na(df[[id_col]]) & nzchar(as.character(df[[id_col]]))
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with missing row labels.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

for (col in value_cols) {
    df[[col]] <- as.numeric(df[[col]])
}

mat <- as.matrix(as.data.frame(df[value_cols], optional = TRUE))
rownames(mat) <- as.character(df[[id_col]])

k_row <- as.integer(config$cutree_rows)
k_col <- as.integer(config$cutree_cols)
if (!is.finite(k_row) || k_row < 1L) {
    k_row <- 1L
}
if (!is.finite(k_col) || k_col < 1L) {
    k_col <- 1L
}
k_row <- min(k_row, nrow(mat))
k_col <- min(k_col, ncol(mat))

row_hc <- stats::as.dendrogram(
    stats::hclust(stats::dist(mat), method = config$clustering_method)
)
col_hc <- stats::as.dendrogram(
    stats::hclust(stats::dist(t(mat)), method = config$clustering_method)
)
row_split <- if (k_row > 1L) k_row else NULL
col_split <- if (k_col > 1L) k_col else NULL

fill_lim <- as.numeric(config$fill_limits)
if (length(fill_lim) != 2L || any(!is.finite(fill_lim))) {
    abs_max <- max(abs(mat), na.rm = TRUE)
    if (!is.finite(abs_max) || abs_max == 0) {
        abs_max <- 1
    }
    fill_lim <- c(-abs_max, abs_max)
}

fill_cols <- palette_colors("Diverging.Tropic")
n_fill <- length(fill_cols)
col_fun <- circlize::colorRamp2(
    seq(fill_lim[[1]], fill_lim[[2]], length.out = n_fill),
    fill_cols
)
na_col <- palette_colors("Qualitative.Safe")[[12]]
border_col <- na_col

left_anno <- NULL
lgd_group <- NULL
row_path <- file.path(dirname(io$input), "rowInfo.tsv")
if (file.exists(row_path)) {
    row_info <- read_table_auto(row_path)
    row_id <- names(row_info)[[1]]
    group_col <- setdiff(names(row_info), row_id)[[1]]
    matched <- match(rownames(mat), as.character(row_info[[row_id]]))
    if (!all(is.na(matched))) {
        groups <- as.character(row_info[[group_col]][matched])
        names(groups) <- rownames(mat)
        levels <- unique(groups[!is.na(groups)])
        named <- config$group_colors
        named <- named[names(named) %in% levels]
        leftover <- setdiff(levels, names(named))
        if (length(leftover)) {
            extra <- palette_colors("Qualitative.Bold", n = length(leftover))
            names(extra) <- leftover
            named <- c(named, extra)
        }
        fill_map <- named[levels]
        left_anno <- ComplexHeatmap::rowAnnotation(
            Group = groups,
            col = list(Group = fill_map),
            show_annotation_name = FALSE,
            show_legend = FALSE,
            simple_anno_size = grid::unit(3.2, "mm")
        )
        lgd_group <- ComplexHeatmap::Legend(
            labels = names(fill_map),
            legend_gp = grid::gpar(fill = unname(fill_map)),
            title = NULL,
            direction = "horizontal",
            nrow = 1,
            grid_width = grid::unit(3.5, "mm"),
            grid_height = grid::unit(3.5, "mm"),
            gap = grid::unit(2, "mm"),
            labels_gp = grid::gpar(fontsize = 8)
        )
    }
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
ht <- ComplexHeatmap::Heatmap(
    mat,
    name = config$labels$fill,
    col = col_fun,
    na_col = na_col,
    cluster_rows = row_hc,
    cluster_columns = col_hc,
    row_split = row_split,
    column_split = col_split,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    row_gap = grid::unit(0.5, "mm"),
    column_gap = grid::unit(0.5, "mm"),
    row_title = NULL,
    column_title = NULL,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "right",
    column_names_rot = 90,
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 9),
    row_dend_width = grid::unit(12, "pt"),
    column_dend_height = grid::unit(10, "pt"),
    left_annotation = left_anno,
    rect_gp = grid::gpar(col = border_col, lwd = 0.35),
    border = FALSE,
    width = ncol(mat) * grid::unit(10, "pt"),
    height = nrow(mat) * grid::unit(10, "pt"),
    show_heatmap_legend = FALSE
)

lgd_fill <- ComplexHeatmap::Legend(
    col_fun = col_fun,
    title = config$labels$fill,
    at = c(fill_lim[[1]], 0, fill_lim[[2]]),
    direction = "horizontal",
    legend_width = grid::unit(3.4, "cm"),
    title_position = "topcenter"
)
lgd <- if (is.null(lgd_group)) {
    lgd_fill
} else {
    packed <- ComplexHeatmap::packLegend(
        lgd_group,
        lgd_fill,
        direction = "vertical",
        gap = grid::unit(3, "mm")
    )
    parent_w <- packed@grob$vp$width
    for (nm in names(packed@grob$children)) {
        child_vp <- packed@grob$children[[nm]]$vp
        packed@grob$children[[nm]]$vp$x <- (parent_w - child_vp$width) * 0.5
    }
    packed
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
width <- 9.4
height <- 8.6
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
tryCatch(
    ComplexHeatmap::draw(
        ht,
        annotation_legend_list = list(lgd),
        annotation_legend_side = "bottom",
        align_annotation_legend = "heatmap_center",
        column_title = config$labels$title,
        column_title_side = "top",
        row_title = config$labels$y,
        row_title_side = "left",
        padding = grid::unit(c(6, 4, 2, 2), "mm"),
        background = NA
    ),
    finally = grDevices::dev.off()
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
