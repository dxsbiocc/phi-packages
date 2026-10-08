#!/usr/bin/env Rscript

# Template-ID: heatmap-oncoprint
#
# Purpose:
#   Draw a cBioPortal-style oncoprint: genes × samples, one cell
#   holding layered alteration glyphs (CNV fill plus a mutation bar),
#   with frequency bars on the top and right.
#
# Inputs:
#   One row per gene-sample-alteration. Default example:
#     - gene, sample, alteration (MUT / AMP / HOMDEL)
#   Several alterations in one cell are several rows, or one cell
#   with MUT;AMP.
#
# Output:
#   A PDF, PNG, or SVG oncoprint heatmap.
#
# Dependencies:
#   ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Extra alteration types are extra rows and extra names in
#   alter_order / bar_alterations, not extra ids. Triangle vs bar
#   glyphs, bar side, pct side, and empty-row dropping stay in
#   CONFIG. An expression matrix beside the oncoprint is a sidecar
#   or a multi-panel figure, not a second id. Do not parse MAF,
#   fetch cBioPortal, or call variants in this script.
#
# Scientific assumptions:
#   Alterations are supplied. Default column/row order is oncoPrint's
#   recurrence sort of those events, a display order, not a test.
#   Percentages are the fraction of samples with any alteration in
#   that gene among the plotted columns.

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
        sample = "sample",
        alteration = "alteration"
    ),
    alter_order = c("HOMDEL", "AMP", "MUT"),
    bar_alterations = "MUT",
    palettes = list(
        alteration = "Qualitative.Safe"
    ),
    labels = list(
        title = "",
        legend = "Alteration",
        alterations = list(
            HOMDEL = "Deep deletion",
            AMP = "Amplification",
            MUT = "Mutation"
        )
    ),
    sort_by_frequency = TRUE,
    remove_empty_rows = TRUE,
    remove_empty_columns = TRUE,
    show_column_names = TRUE,
    show_pct = TRUE,
    row_names_side = "left",
    pct_side = "right",
    show_top_bar = TRUE,
    show_right_bar = TRUE,
    size = list(
        width = 9.4,
        height = 5.8
    )
)

load_packages(c("ComplexHeatmap", "grid", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
g_col <- config$columns$gene
s_col <- config$columns$sample
a_col <- config$columns$alteration
df[[g_col]] <- as.character(df[[g_col]])
df[[s_col]] <- as.character(df[[s_col]])
df[[a_col]] <- as.character(df[[a_col]])
ok <- !is.na(df[[g_col]]) & nzchar(df[[g_col]]) &
    !is.na(df[[s_col]]) & nzchar(df[[s_col]]) &
    !is.na(df[[a_col]]) & nzchar(trimws(df[[a_col]]))
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with missing gene, sample, or alteration.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No alterations left to plot.", call. = FALSE)
}
df[[a_col]] <- trimws(df[[a_col]])
parts <- strsplit(df[[a_col]], "\\s*[;:,|]\\s*")
df <- data.frame(
    gene = rep(df[[g_col]], lengths(parts)),
    sample = rep(df[[s_col]], lengths(parts)),
    alteration = trimws(unlist(parts, use.names = FALSE)),
    stringsAsFactors = FALSE
)
df <- df[nzchar(df$alteration), , drop = FALSE]
df <- unique(df)
if (!nrow(df)) {
    stop("No alterations left after splitting type strings.", call. = FALSE)
}

genes <- unique(df$gene)
samples <- unique(df$sample)
type_order <- unique(c(
    intersect(as.character(config$alter_order), df$alteration),
    setdiff(unique(df$alteration), config$alter_order)
))
if (!length(type_order)) {
    stop("No alteration types to draw.", call. = FALSE)
}

mat <- matrix(
    "",
    nrow = length(genes),
    ncol = length(samples),
    dimnames = list(genes, samples)
)
cell <- split(df$alteration, paste(df$gene, df$sample, sep = "\r"))
for (nm in names(cell)) {
    xy <- strsplit(nm, "\r", fixed = TRUE)[[1]]
    types <- unique(cell[[nm]])
    types <- types[types %in% type_order]
    types <- type_order[type_order %in% types]
    mat[xy[[1]], xy[[2]]] <- paste(types, collapse = ";")
}

pal_all <- palette_colors(config$palettes$alteration)
if (identical(config$palettes$alteration, "Qualitative.Safe") &&
        length(pal_all) >= 4L) {
    # Skip gold so default HOMDEL/AMP/MUT are cyan, rose, green.
    pick <- pal_all[c(1, 2, 4, setdiff(seq_along(pal_all), c(1, 2, 4)))]
} else {
    pick <- pal_all
}
alter_col <- stats::setNames(
    expand_palette(pick, length(type_order)),
    type_order
)
bg <- grDevices::adjustcolor(
    palette_colors("Qualitative.Safe")[[12]],
    alpha.f = 0.2
)
bar_types <- intersect(as.character(config$bar_alterations), type_order)
full_types <- setdiff(type_order, bar_types)

make_rect <- function(fill, height) {
    force(fill)
    force(height)
    function(x, y, w, h) {
        grid::grid.rect(
            x, y, w * 0.9, h * height,
            gp = grid::gpar(fill = fill, col = NA)
        )
    }
}
alter_fun <- list(
    background = make_rect(bg, 0.9)
)
for (tp in full_types) {
    alter_fun[[tp]] <- make_rect(unname(alter_col[[tp]]), 0.9)
}
for (tp in bar_types) {
    alter_fun[[tp]] <- make_rect(unname(alter_col[[tp]]), 0.33)
}

legend_labels <- vapply(type_order, function(tp) {
    lab <- config$labels$alterations[[tp]]
    if (is.null(lab) || !nzchar(as.character(lab))) {
        tp
    } else {
        as.character(lab)
    }
}, character(1), USE.NAMES = TRUE)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
top_anno <- if (isTRUE(config$show_top_bar)) {
    ComplexHeatmap::HeatmapAnnotation(
        cbar = ComplexHeatmap::anno_oncoprint_barplot(
            height = grid::unit(1.4, "cm")
        ),
        show_annotation_name = FALSE
    )
} else {
    NULL
}
right_anno <- if (isTRUE(config$show_right_bar)) {
    ComplexHeatmap::rowAnnotation(
        rbar = ComplexHeatmap::anno_oncoprint_barplot(
            width = grid::unit(1.6, "cm")
        ),
        show_annotation_name = FALSE
    )
} else {
    NULL
}

ht_args <- list(
    mat = mat,
    alter_fun = alter_fun,
    col = alter_col,
    remove_empty_columns = isTRUE(config$remove_empty_columns),
    remove_empty_rows = isTRUE(config$remove_empty_rows),
    show_column_names = isTRUE(config$show_column_names),
    show_pct = isTRUE(config$show_pct),
    row_names_side = config$row_names_side,
    pct_side = config$pct_side,
    top_annotation = top_anno,
    right_annotation = right_anno,
    column_title = if (nzchar(config$labels$title)) {
        config$labels$title
    } else {
        NULL
    },
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 8),
    pct_gp = grid::gpar(fontsize = 8),
    heatmap_legend_param = list(
        title = config$labels$legend,
        at = type_order,
        labels = unname(legend_labels[type_order]),
        nrow = 1,
        title_position = "leftcenter"
    ),
    alter_fun_is_vectorized = FALSE
)
if (!isTRUE(config$sort_by_frequency)) {
    ht_args$row_order <- seq_len(nrow(mat))
    ht_args$column_order <- seq_len(ncol(mat))
}

ht <- do.call(ComplexHeatmap::oncoPrint, ht_args)

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
    ComplexHeatmap::draw(
        ht,
        heatmap_legend_side = "bottom",
        padding = grid::unit(c(2, 2, 2, 2), "mm"),
        background = NA
    ),
    finally = grDevices::dev.off()
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
