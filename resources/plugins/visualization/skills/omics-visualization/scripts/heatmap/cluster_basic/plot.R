#!/usr/bin/env Rscript

# Template-ID: heatmap-cluster-basic
#
# Purpose:
#   Draw a clustered heatmap of a feature-by-sample numeric matrix,
#   with optional row and column annotations.
#
# Inputs:
#   A table with one feature-id column and numeric sample columns.
#   Default example:
#     - sample: feature (gene) identifier
#     - remaining columns: sample expression values
#   Optional tables beside the input (same directory):
#     - colInfo.tsv: sample, group
#     - rowInfo.tsv: gene, group
#
# Output:
#   A PDF, PNG, or SVG clustered heatmap.
#
# Dependencies:
#   ComplexHeatmap, circlize, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Place colInfo.tsv and rowInfo.tsv next to the input if annotations
#   are needed; the CLI stays one input plus one output.
#   Edit DATA PREPARATION to change scaling or clustering.
#   Edit PLOT only when the heatmap geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one feature and each numeric column is one sample.
#   Rows are z-scored (center and scale) for display; this is not an
#   upstream differential-expression test.
#   Dendrograms are hierarchical clustering of the scaled matrix.

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
        id = "sample"
    ),
    labels = list(
        title = "",
        x = "",
        y = "Gene set",
        fill = "z-score"
    )
)

load_packages(c("ComplexHeatmap", "circlize", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

id_col <- config$columns$id
sample_cols <- setdiff(names(df), id_col)
if (length(sample_cols) < 2L) {
    stop("Heatmap input needs at least two numeric sample columns.",
        call. = FALSE
    )
}

for (col in sample_cols) {
    df[[col]] <- as.numeric(df[[col]])
}

mat <- as.matrix(as.data.frame(df[sample_cols], optional = TRUE))
rownames(mat) <- as.character(df[[id_col]])
mat <- t(scale(t(mat), center = TRUE, scale = TRUE))

annotation_colors <- function(anno_df, palette) {
    out <- list()
    for (nm in names(anno_df)) {
        levels <- unique(as.character(anno_df[[nm]]))
        levels <- levels[!is.na(levels)]
        n <- length(levels)
        if (!n) next
        cols <- if (n <= length(palette)) {
            palette[seq_len(n)]
        } else {
            grDevices::colorRampPalette(palette)(n)
        }
        out[[nm]] <- stats::setNames(cols, levels)
    }
    out
}

input_dir <- dirname(io$input)
col_path <- file.path(input_dir, "colInfo.tsv")
row_path <- file.path(input_dir, "rowInfo.tsv")

top_anno <- NULL
if (file.exists(col_path)) {
    col_info <- read_table_auto(col_path)
    col_id <- names(col_info)[[1]]
    col_vars <- setdiff(names(col_info), col_id)
    matched <- match(colnames(mat), col_info[[col_id]])
    if (length(col_vars) && !all(is.na(matched))) {
        anno_df <- as.data.frame(
            col_info[matched, col_vars, drop = FALSE],
            optional = TRUE
        )
        rownames(anno_df) <- colnames(mat)
        top_anno <- HeatmapAnnotation(
            df = anno_df,
            col = annotation_colors(anno_df, c("#7fbf7b", "#af8dc3"))
        )
    }
}

left_anno <- NULL
if (file.exists(row_path)) {
    row_info <- read_table_auto(row_path)
    row_id <- names(row_info)[[1]]
    row_vars <- setdiff(names(row_info), row_id)
    matched <- match(rownames(mat), row_info[[row_id]])
    if (length(row_vars) && !all(is.na(matched))) {
        anno_df <- as.data.frame(
            row_info[matched, row_vars, drop = FALSE],
            optional = TRUE
        )
        rownames(anno_df) <- rownames(mat)
        left_anno <- rowAnnotation(
            df = anno_df,
            col = annotation_colors(anno_df, c("#8c510a", "#01665e"))
        )
    }
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
col_fun <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#8c510a", "white", "#01665e")
)

ht <- Heatmap(
    mat,
    name = config$labels$fill,
    col = col_fun,
    na_col = "transparent",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    top_annotation = top_anno,
    left_annotation = left_anno,
    row_title = config$labels$y,
    column_title = config$labels$title
)

grob <- grid::grid.grabExpr(draw(ht, heatmap_legend_side = "right"))

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(grob, io$output)
