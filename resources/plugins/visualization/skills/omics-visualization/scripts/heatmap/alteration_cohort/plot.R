#!/usr/bin/env Rscript

# Template-ID: heatmap-alteration-cohort
#
# Purpose:
#   Draw a gene × cancer-cohort alteration-frequency heatmap: each
#   cell holds a supplied count coloured by alteration class (with
#   intensity scaled by the value), stacked frequency bars on the
#   top and right, a gene-role strip on the left, and cohorts grouped
#   by organ system along the bottom. An optional Total row sums
#   each column.
#
# Inputs:
#   Main table (one row per non-empty gene × cohort cell):
#     - gene, cohort, alteration, value
#   Sidecars next to the input (or absolute paths):
#     - genes.tsv: gene, role (Oncogene / Tumor Suppressor, …)
#     - cohorts.tsv: cohort, group (organ system; file order = column order)
#
# Output:
#   A PDF, PNG, or SVG ComplexHeatmap figure.
#
# Dependencies:
#   ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Extra alteration classes are extra names in alter_order and
#   palettes$alteration (one "Family.Name:index" per class). Extra
#   gene roles are extra names in role_colors. Gene × sample layered
#   glyphs stay heatmap-oncoprint. Do not call variants, fetch
#   cBioPortal, or run enrichment in this script.
#
# Scientific assumptions:
#   Values, alteration class per cell, gene roles and cohort groups
#   are supplied. Intensity is a display scale of the supplied
#   values; bar heights are sums of those values. No test is run.

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
        gene = "gene",
        cohort = "cohort",
        alteration = "alteration",
        value = "value"
    ),
    genes_file = "genes.tsv",
    cohorts_file = "cohorts.tsv",
    gene_columns = list(gene = "gene", role = "role"),
    cohort_columns = list(cohort = "cohort", group = "group"),
    alter_order = c(
        "Amp", "Del", "Fusion", "Mut", "Methylation", "Multiple"
    ),
    # One catalog colour per alteration class (Family.Name:index).
    palettes = list(
        alteration = c(
            Amp = "Qualitative.Set1 / 9:1",
            Del = "Qualitative.Set1 / 9:2",
            Fusion = "Qualitative.Set1 / 9:4",
            Mut = "Qualitative.Set1 / 9:3",
            Methylation = "Qualitative.Bold:4",
            Multiple = "Qualitative.Set1 / 9:5"
        )
    ),
    role_colors = list(
        "Oncogene" = "Qualitative.Set1 / 9:1",
        "Tumor Suppressor" = "Qualitative.Set1 / 9:2"
    ),
    group_palette = "Qualitative.Safe",
    labels = list(
        title = "",
        top_bar = "Alteration frequencies",
        role = "Role",
        legend = "Alteration"
    ),
    show_top_bar = TRUE,
    show_right_bar = TRUE,
    show_role = TRUE,
    show_total = TRUE,
    show_numbers = TRUE,
    empty_fill = "#FFFFFF",
    cell_border = "grey85",
    intensity_floor = 0.18,  # lightest shade for the smallest positive value
    top_bar_height = 2.2,
    right_bar_width = 2.4,
    role_width = 0.35,
    group_height = 0.35,
    fontsize = list(row = 9, column = 8, cell = 7, legend = 9),
    size = list(width = 12.5, height = 8.2)
)

load_packages(c("ComplexHeatmap", "circlize", "grid", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

sidecar_path <- function(name) {
    if (is.null(name) || !nzchar(name)) return(NULL)
    path <- if (file.exists(name)) name else file.path(dirname(io$input), name)
    if (file.exists(path)) path else NULL
}

blend_col <- function(col, t, bg = config$empty_fill) {
    t <- max(0, min(1, as.numeric(t)))
    rgb <- grDevices::col2rgb(col) / 255
    bgv <- grDevices::col2rgb(bg) / 255
    grDevices::rgb(
        (1 - t) * bgv[1] + t * rgb[1],
        (1 - t) * bgv[2] + t * rgb[2],
        (1 - t) * bgv[3] + t * rgb[3]
    )
}

df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
g_col <- config$columns$gene
c_col <- config$columns$cohort
a_col <- config$columns$alteration
v_col <- config$columns$value
df[[g_col]] <- as.character(df[[g_col]])
df[[c_col]] <- as.character(df[[c_col]])
df[[a_col]] <- as.character(df[[a_col]])
df[[v_col]] <- as.numeric(df[[v_col]])
ok <- !is.na(df[[g_col]]) & nzchar(df[[g_col]]) &
    !is.na(df[[c_col]]) & nzchar(df[[c_col]]) &
    !is.na(df[[a_col]]) & nzchar(df[[a_col]]) &
    is.finite(df[[v_col]]) & df[[v_col]] > 0
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with missing fields or non-positive value.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) stop("No cells left to plot.", call. = FALSE)

genes_path <- sidecar_path(config$genes_file)
cohorts_path <- sidecar_path(config$cohorts_file)
if (is.null(genes_path)) stop("Missing genes sidecar: ", config$genes_file, call. = FALSE)
if (is.null(cohorts_path)) {
    stop("Missing cohorts sidecar: ", config$cohorts_file, call. = FALSE)
}
genes_df <- as.data.frame(read_table_auto(genes_path), stringsAsFactors = FALSE)
cohorts_df <- as.data.frame(read_table_auto(cohorts_path), stringsAsFactors = FALSE)
require_columns(genes_df, config$gene_columns)
require_columns(cohorts_df, config$cohort_columns)
genes_df[[config$gene_columns$gene]] <- as.character(
    genes_df[[config$gene_columns$gene]]
)
genes_df[[config$gene_columns$role]] <- as.character(
    genes_df[[config$gene_columns$role]]
)
cohorts_df[[config$cohort_columns$cohort]] <- as.character(
    cohorts_df[[config$cohort_columns$cohort]]
)
cohorts_df[[config$cohort_columns$group]] <- as.character(
    cohorts_df[[config$cohort_columns$group]]
)

gene_order <- unique(genes_df[[config$gene_columns$gene]])
cohort_order <- unique(cohorts_df[[config$cohort_columns$cohort]])
role_map <- stats::setNames(
    genes_df[[config$gene_columns$role]],
    genes_df[[config$gene_columns$gene]]
)
group_map <- stats::setNames(
    cohorts_df[[config$cohort_columns$group]],
    cohorts_df[[config$cohort_columns$cohort]]
)

# Keep only genes / cohorts listed in the sidecars; append any extras at end.
extra_g <- setdiff(unique(df[[g_col]]), gene_order)
extra_c <- setdiff(unique(df[[c_col]]), cohort_order)
if (length(extra_g)) {
    message("Appending ", length(extra_g), " gene(s) missing from genes.tsv.")
    gene_order <- c(gene_order, extra_g)
    role_map[extra_g] <- NA_character_
}
if (length(extra_c)) {
    message("Appending ", length(extra_c), " cohort(s) missing from cohorts.tsv.")
    cohort_order <- c(cohort_order, extra_c)
    group_map[extra_c] <- "Other"
}
df <- df[df[[g_col]] %in% gene_order & df[[c_col]] %in% cohort_order, , drop = FALSE]

type_order <- unique(c(
    intersect(as.character(config$alter_order), df[[a_col]]),
    setdiff(unique(df[[a_col]]), config$alter_order)
))
if (!length(type_order)) stop("No alteration types to draw.", call. = FALSE)

alter_col <- vapply(type_order, function(tp) {
    ref <- config$palettes$alteration[[tp]]
    if (is.null(ref)) {
        stop("No palette entry for alteration: ", tp, call. = FALSE)
    }
    pal_ref(ref)
}, character(1), USE.NAMES = TRUE)

val_mat <- matrix(
    0, nrow = length(gene_order), ncol = length(cohort_order),
    dimnames = list(gene_order, cohort_order)
)
type_mat <- matrix(
    NA_character_, nrow = length(gene_order), ncol = length(cohort_order),
    dimnames = list(gene_order, cohort_order)
)
for (i in seq_len(nrow(df))) {
    g <- df[[g_col]][i]
    c <- df[[c_col]][i]
    val_mat[g, c] <- df[[v_col]][i]
    type_mat[g, c] <- df[[a_col]][i]
}

vmax <- max(val_mat, na.rm = TRUE)
if (!is.finite(vmax) || vmax <= 0) stop("All values are zero.", call. = FALSE)

# Stacked bar matrices: types × cohorts (top) and genes × types (right).
top_mat <- matrix(
    0, nrow = length(type_order), ncol = length(cohort_order),
    dimnames = list(type_order, cohort_order)
)
right_mat <- matrix(
    0, nrow = length(gene_order), ncol = length(type_order),
    dimnames = list(gene_order, type_order)
)
for (i in seq_len(nrow(df))) {
    g <- df[[g_col]][i]
    c <- df[[c_col]][i]
    a <- df[[a_col]][i]
    v <- df[[v_col]][i]
    top_mat[a, c] <- top_mat[a, c] + v
    right_mat[g, a] <- right_mat[g, a] + v
}

roles <- unname(role_map[gene_order])
groups <- unname(group_map[cohort_order])
group_levels <- unique(groups)
group_cols <- expand_palette(
    palette_colors(config$group_palette), length(group_levels)
)
names(group_cols) <- group_levels

role_levels <- unique(na.omit(roles))
role_fill <- vapply(role_levels, function(r) {
    ref <- config$role_colors[[r]]
    if (is.null(ref)) palette_colors("Qualitative.Safe")[[12]] else pal_ref(ref)
}, character(1), USE.NAMES = TRUE)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
cell_fun <- function(j, i, x, y, w, h, fill) {
    tp <- type_mat[i, j]
    v <- val_mat[i, j]
    if (!is.finite(v) || v <= 0 || is.na(tp)) {
        grid::grid.rect(
            x, y, w, h,
            gp = grid::gpar(fill = config$empty_fill, col = config$cell_border)
        )
        return(invisible())
    }
    t <- config$intensity_floor +
        (1 - config$intensity_floor) * (v / vmax)
    fill_col <- blend_col(alter_col[[tp]], t)
    grid::grid.rect(
        x, y, w, h,
        gp = grid::gpar(fill = fill_col, col = config$cell_border)
    )
    if (isTRUE(config$show_numbers)) {
        grid::grid.text(
            as.character(round(v)), x, y,
            gp = grid::gpar(
                fontsize = config$fontsize$cell,
                col = if (t > 0.72) "white" else "grey15"
            )
        )
    }
}

top_anno <- NULL
if (isTRUE(config$show_top_bar)) {
    top_anno <- ComplexHeatmap::HeatmapAnnotation(
        `Alteration frequencies` = ComplexHeatmap::anno_barplot(
            t(top_mat),
            gp = grid::gpar(fill = alter_col[type_order], col = NA),
            border = FALSE,
            bar_width = 0.7,
            height = grid::unit(config$top_bar_height, "cm")
        ),
        annotation_name_side = "left",
        annotation_name_gp = grid::gpar(fontsize = 8),
        show_legend = FALSE
    )
}

bottom_anno <- ComplexHeatmap::HeatmapAnnotation(
    group = groups,
    col = list(group = group_cols),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontsize = 8),
    annotation_legend_param = list(
        group = list(title = "Organ system", nrow = 2)
    ),
    simple_anno_size = grid::unit(config$group_height, "cm"),
    show_legend = TRUE
)

left_anno <- NULL
if (isTRUE(config$show_role)) {
    left_anno <- ComplexHeatmap::rowAnnotation(
        Role = roles,
        col = list(Role = role_fill),
        annotation_name_side = "top",
        annotation_name_gp = grid::gpar(fontsize = 8),
        simple_anno_size = grid::unit(config$role_width, "cm"),
        annotation_legend_param = list(
            Role = list(title = config$labels$role, nrow = 1)
        )
    )
}

right_anno <- NULL
if (isTRUE(config$show_right_bar)) {
    right_anno <- ComplexHeatmap::rowAnnotation(
        freq = ComplexHeatmap::anno_barplot(
            right_mat,
            gp = grid::gpar(fill = alter_col[type_order], col = NA),
            border = FALSE,
            bar_width = 0.7,
            width = grid::unit(config$right_bar_width, "cm"),
            axis_param = list(side = "top")
        ),
        show_annotation_name = FALSE
    )
}

# Dummy continuous col so Heatmap accepts a numeric matrix; cells are painted
# entirely by cell_fun.
ht <- ComplexHeatmap::Heatmap(
    val_mat,
    name = "value",
    col = circlize::colorRamp2(c(0, vmax), c(config$empty_fill, "grey50")),
    cell_fun = cell_fun,
    rect_gp = grid::gpar(type = "none"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = !isTRUE(config$show_total),
    row_names_side = "left",
    column_names_rot = 90,
    column_names_centered = FALSE,
    row_names_gp = grid::gpar(fontsize = config$fontsize$row),
    column_names_gp = grid::gpar(
        fontsize = config$fontsize$column,
        col = group_cols[groups]
    ),
    column_split = factor(groups, levels = group_levels),
    cluster_column_slices = FALSE,
    column_title = if (isTRUE(config$show_total)) NULL else character(0),
    column_title_side = "bottom",
    column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    column_gap = grid::unit(0.8, "mm"),
    top_annotation = top_anno,
    bottom_annotation = if (isTRUE(config$show_total)) NULL else bottom_anno,
    left_annotation = left_anno,
    right_annotation = right_anno,
    show_heatmap_legend = FALSE,
    border = FALSE
)

# Total row as a second heatmap under the main one.
ht_total <- NULL
if (isTRUE(config$show_total)) {
    total_vec <- matrix(
        colSums(val_mat), nrow = 1,
        dimnames = list("Total", cohort_order)
    )
    total_type <- matrix(
        NA_character_, nrow = 1, ncol = length(cohort_order),
        dimnames = list("Total", cohort_order)
    )
    # Colour Total by the dominant alteration in that column.
    for (j in seq_along(cohort_order)) {
        col_vals <- top_mat[, j]
        if (sum(col_vals) > 0) {
            total_type[1, j] <- names(col_vals)[which.max(col_vals)]
        }
    }
    cell_fun_total <- function(j, i, x, y, w, h, fill) {
        tp <- total_type[i, j]
        v <- total_vec[i, j]
        if (!is.finite(v) || v <= 0 || is.na(tp)) {
            grid::grid.rect(
                x, y, w, h,
                gp = grid::gpar(fill = config$empty_fill, col = config$cell_border)
            )
            return(invisible())
        }
        t <- config$intensity_floor +
            (1 - config$intensity_floor) * (v / max(total_vec))
        fill_col <- blend_col(alter_col[[tp]], t)
        grid::grid.rect(
            x, y, w, h,
            gp = grid::gpar(fill = fill_col, col = config$cell_border)
        )
        if (isTRUE(config$show_numbers)) {
            grid::grid.text(
                as.character(round(v)), x, y,
                gp = grid::gpar(
                    fontsize = config$fontsize$cell,
                    fontface = "bold",
                    col = if (t > 0.72) "white" else "grey15"
                )
            )
        }
    }
    right_total <- NULL
    if (isTRUE(config$show_right_bar)) {
        right_total <- ComplexHeatmap::rowAnnotation(
            freq = ComplexHeatmap::anno_empty(
                width = grid::unit(config$right_bar_width, "cm"),
                border = FALSE
            ),
            show_annotation_name = FALSE
        )
    }
    left_total <- NULL
    if (isTRUE(config$show_role)) {
        left_total <- ComplexHeatmap::rowAnnotation(
            Role = ComplexHeatmap::anno_empty(
                width = grid::unit(config$role_width, "cm"),
                border = FALSE
            ),
            show_annotation_name = FALSE
        )
    }
    ht_total <- ComplexHeatmap::Heatmap(
        total_vec,
        name = "total",
        col = circlize::colorRamp2(c(0, max(total_vec)), c(config$empty_fill, "grey50")),
        cell_fun = cell_fun_total,
        rect_gp = grid::gpar(type = "none"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = config$fontsize$row, fontface = "bold"),
        column_names_rot = 90,
        column_names_gp = grid::gpar(
            fontsize = config$fontsize$column,
            col = group_cols[groups]
        ),
        column_split = factor(groups, levels = group_levels),
        cluster_column_slices = FALSE,
        column_title_side = "bottom",
        column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
        column_gap = grid::unit(0.8, "mm"),
        bottom_annotation = bottom_anno,
        left_annotation = left_total,
        right_annotation = right_total,
        show_heatmap_legend = FALSE,
        border = FALSE,
        height = grid::unit(6, "mm")
    )
}

alter_legend <- ComplexHeatmap::Legend(
    title = config$labels$legend,
    at = type_order,
    legend_gp = grid::gpar(fill = unname(alter_col[type_order])),
    nrow = 1,
    title_position = "leftcenter",
    labels_gp = grid::gpar(fontsize = config$fontsize$legend),
    title_gp = grid::gpar(fontsize = config$fontsize$legend, fontface = "bold")
)

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
width <- as.numeric(config$size$width)
height <- as.numeric(config$size$height)
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
    {
        ht_list <- if (is.null(ht_total)) ht else ht %v% ht_total
        ComplexHeatmap::draw(
            ht_list,
            annotation_legend_list = list(alter_legend),
            annotation_legend_side = "bottom",
            heatmap_legend_side = "bottom",
            merge_legend = TRUE,
            padding = grid::unit(c(4, 2, 2, 6), "mm"),
            background = NA,
            column_title = if (nzchar(config$labels$title)) {
                config$labels$title
            } else {
                NULL
            }
        )
    },
    finally = grDevices::dev.off()
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
