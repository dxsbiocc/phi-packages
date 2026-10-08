#!/usr/bin/env Rscript

# Template-ID: heatmap-enrichment-zoom
#
# Purpose:
#   Draw an unclustered DE expression heatmap split by sample group and
#   gene direction, with optional TF / surface tracks and enrichment
#   bar panels (anno_zoom) aligned to the Up / Down row slices.
#   Extra databases are extra columns in enrichment.tsv, not extra ids.
#   Right-hand zoom panels are ggplot2 bar charts (geom_col) printed into
#   the anno_zoom viewport. Do not replace them with grid.text lists.
#
# Inputs:
#   A table with one gene-id column and numeric sample columns.
#   Default example:
#     - gene: feature identifier
#     - remaining columns: sample expression values
#   Optional tables beside the input (same directory):
#     - colInfo.tsv: sample, group
#     - rowInfo.tsv: gene, direct, optional TF, SP, mark
#     - enrichment.tsv: database, Change, Description, p.adjust
#   Change must match row-split levels (Up / Down). p.adjust is a
#   supplied enrichment result; this script does not run enricher.
#
# Output:
#   A PDF, PNG, or SVG composite heatmap with zoomed enrichment bars.
#
# Dependencies:
#   ComplexHeatmap, circlize, ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names, n_terms, size).
#   Place sidecars next to the input; the CLI stays one input plus one
#   output. Edit DATA PREPARATION to change scaling or split order.
#   Edit PLOT only when the heatmap or zoom geometry must change.
#   Ranked Cartesian enrichment alone is bar-enrichment-*. A clustered
#   matrix without zoom panels is heatmap-cluster-basic. Do not add an
#   id per database (GO / KEGG / WikiPathway). Zoom panels are ggplot2
#   bar charts, not a second template.
#
# Scientific assumptions:
#   The matrix is already a DE gene set; this script does not run
#   limma or an equivalent test. Rows are z-scored for display.
#   Enrichment terms and p.adjust are supplied. -log10(p.adjust) is a
#   display transform. Clustering is off so row/column order is the
#   supplied (split) order, not a statistical grouping.

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
        id = "gene"
    ),
    sidecars = list(
        col = "colInfo.tsv",
        row = "rowInfo.tsv",
        enrichment = "enrichment.tsv"
    ),
    row = list(
        split = "direct",
        tf = "TF",
        surface = "SP",
        mark = "mark"
    ),
    col = list(
        split = "group"
    ),
    split_levels = list(
        row = c("Down", "Up"),
        column = c("Normal", "Tumor")
    ),
    n_terms = 4L,
    zoom = list(
        height_cm = 4,
        width_cm = 8
    ),
    labels = list(
        fill = "zscore",
        direct = "Direct",
        tf = "TF",
        surface = "SP",
        p_up = "-log10(Pvalue)-Up",
        p_down = "-log10(Pvalue)-Down"
    ),
    size = list(
        width = 18,
        height = 6
    )
)

load_packages(c("ComplexHeatmap", "circlize", "ggplot2", "readr"))

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

input_dir <- dirname(io$input)
col_path <- file.path(input_dir, config$sidecars$col)
row_path <- file.path(input_dir, config$sidecars$row)
enr_path <- file.path(input_dir, config$sidecars$enrichment)

if (!file.exists(col_path) || !file.exists(row_path) || !file.exists(enr_path)) {
    stop(
        paste(
            "This template needs colInfo.tsv, rowInfo.tsv, and",
            "enrichment.tsv beside the input matrix."
        ),
        call. = FALSE
    )
}

col_info <- read_table_auto(col_path)
row_info <- read_table_auto(row_path)
enr <- read_table_auto(enr_path)
require_columns(enr, list(
    database = "database",
    Change = "Change",
    Description = "Description",
    p.adjust = "p.adjust"
))

col_id <- names(col_info)[[1]]
row_id <- names(row_info)[[1]]
col_matched <- match(colnames(mat), col_info[[col_id]])
row_matched <- match(rownames(mat), row_info[[row_id]])
if (anyNA(col_matched) || anyNA(row_matched)) {
    stop("Every sample and gene in the matrix must appear in the sidecars.",
        call. = FALSE
    )
}

group <- as.character(col_info[[config$col$split]][col_matched])
direct <- as.character(row_info[[config$row$split]][row_matched])
group <- factor(group, levels = config$split_levels$column)
direct <- factor(direct, levels = config$split_levels$row)
if (anyNA(group) || anyNA(direct)) {
    stop("split values must be among config$split_levels.", call. = FALSE)
}

ord <- order(direct, rownames(mat))
mat <- mat[ord, , drop = FALSE]
direct <- droplevels(direct[ord])
row_info <- row_info[row_matched[ord], , drop = FALSE]

col_ord <- order(group, colnames(mat))
mat <- mat[, col_ord, drop = FALSE]
group <- droplevels(group[col_ord])

as_flag <- function(x) {
    x <- as.character(x)
    toupper(x) %in% c("TRUE", "T", "1", "YES")
}

mark <- if (config$row$mark %in% names(row_info)) {
    as_flag(row_info[[config$row$mark]])
} else {
    rep(FALSE, nrow(mat))
}

tf_vec <- if (config$row$tf %in% names(row_info)) {
    x <- as.character(row_info[[config$row$tf]])
    x[x %in% c("", "NA")] <- NA_character_
    x
} else {
    NULL
}

sp_vec <- if (config$row$surface %in% names(row_info)) {
    x <- as.character(row_info[[config$row$surface]])
    x[x %in% c("", "NA")] <- NA_character_
    x
} else {
    NULL
}

enr$database <- as.character(enr$database)
enr$Change <- as.character(enr$Change)
enr$Description <- as.character(enr$Description)
enr$p.adjust <- as.numeric(enr$p.adjust)
enr <- enr[is.finite(enr$p.adjust) & enr$p.adjust > 0, , drop = FALSE]

n_keep <- as.integer(config$n_terms)
enr <- do.call(rbind, lapply(split(enr, list(enr$database, enr$Change), drop = TRUE), function(part) {
    part[seq_len(min(n_keep, nrow(part))), , drop = FALSE]
}))
rownames(enr) <- NULL

dbs <- unique(enr$database)
if (!length(dbs)) {
    stop("enrichment.tsv has no usable rows after filtering p.adjust.",
        call. = FALSE
    )
}

neg_log10 <- function(p) {
    p <- as.numeric(p)
    out <- -log10(p)
    out[!is.finite(out)] <- NA_real_
    out
}

maxP <- max(neg_log10(enr$p.adjust), na.rm = TRUE)
if (!is.finite(maxP) || maxP <= 0) {
    maxP <- 1
}

# Colours match the supplied figure. Direct / SP / TF cofactor hex also
# sit in colors.json (Brand.Asana, Brand.sprite, Brand.ns). Heatmap
# endpoints and TF yellow are the figure's hex.
regulate_cols <- c("#01cd74", "#ff4e00")
direct_cols <- c(Down = "#3be8b0", Up = "#fc636b")
tf_cols <- c(TF = "#ffc917", `TF cofactor` = "#003082")
sp_cols <- c(Surface = "#6bc048")

col_fun <- circlize::colorRamp2(
    c(-2, 0, 2),
    c(regulate_cols[[1]], "white", regulate_cols[[2]])
)

ramp_pair <- function(high) {
    pal <- grDevices::colorRampPalette(c("white", high))(6)
    pal[c(2L, 6L)]
}

gradient_up <- circlize::colorRamp2(c(0, maxP), ramp_pair(direct_cols[["Up"]]))
gradient_down <- circlize::colorRamp2(
    c(0, maxP),
    ramp_pair(direct_cols[["Down"]])
)

# Same as the source ggplot_zoom: bars + italic labels, no axes.
theme_zoom <- ggplot2::theme_void() +
    ggplot2::theme(
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank(),
        plot.margin = grid::unit(c(0, 0, 0, 0), "cm")
    )

ggplot_zoom <- function(slice_name, db_tbl, db_name) {
    sub <- db_tbl[db_tbl$Change == slice_name, , drop = FALSE]
    grid::pushViewport(grid::viewport())
    if (identical(slice_name, config$split_levels$row[[1]])) {
        grid::grid.text(
            db_name,
            y = grid::unit(1, "npc") + grid::unit(4, "mm"),
            gp = grid::gpar(fontsize = 11, fontface = "bold")
        )
    }
    grid::grid.rect()
    if (!nrow(sub)) {
        grid::popViewport()
        return(invisible())
    }
    sub$score <- neg_log10(sub$p.adjust)
    sub$Description <- factor(
        sub$Description,
        levels = rev(unique(sub$Description))
    )
    fill_fun <- if (identical(slice_name, "Up")) gradient_up else gradient_down
    p <- ggplot2::ggplot(sub, ggplot2::aes(x = score, y = Description)) +
        ggplot2::geom_col(
            ggplot2::aes(fill = score),
            show.legend = FALSE
        ) +
        ggplot2::geom_text(
            ggplot2::aes(x = 0.1, label = Description),
            hjust = 0,
            fontface = "italic",
            show.legend = FALSE
        ) +
        ggplot2::scale_x_continuous(
            expand = c(0.01, 0.01),
            limits = c(0, maxP)
        ) +
        ggplot2::scale_fill_gradientn(
            colours = fill_fun(seq(0, maxP, length.out = 32L))
        ) +
        theme_zoom
    print(p, newpage = FALSE)
    grid::popViewport()
    invisible()
}

make_zoom <- function(db_name, db_tbl) {
    force(db_name)
    force(db_tbl)
    function(index, name) {
        ggplot_zoom(name, db_tbl, db_name)
    }
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
left_args <- list()
if (any(mark)) {
    left_args$gene <- ComplexHeatmap::anno_mark(
        at = which(mark),
        labels = rownames(mat)[mark],
        which = "row",
        side = "left",
        labels_gp = grid::gpar(fontsize = 8)
    )
}
left_args$cluster1 <- ComplexHeatmap::anno_block(
    gp = grid::gpar(fill = unname(direct_cols[levels(direct)])),
    which = "row",
    width = grid::unit(3, "mm")
)
left_args$annotation_label <- c(
    if (any(mark)) c(gene = ""),
    cluster1 = ""
)
left_anno <- do.call(ComplexHeatmap::rowAnnotation, left_args)

right_args <- list()
shown <- character()
if (!is.null(tf_vec) && any(!is.na(tf_vec))) {
    used <- intersect(names(tf_cols), unique(tf_vec[!is.na(tf_vec)]))
    right_args$TF <- ComplexHeatmap::anno_simple(
        tf_vec,
        col = tf_cols[used],
        which = "row",
        na_col = "white",
        border = TRUE
    )
    shown <- c(shown, "TF")
}
if (!is.null(sp_vec) && any(!is.na(sp_vec))) {
    used <- intersect(names(sp_cols), unique(sp_vec[!is.na(sp_vec)]))
    right_args$SP <- ComplexHeatmap::anno_simple(
        sp_vec,
        col = sp_cols[used],
        which = "row",
        na_col = "white",
        border = TRUE
    )
    shown <- c(shown, "SP")
}

slice_n <- nlevels(direct)
zoom_size <- grid::unit(rep(config$zoom$height_cm, slice_n), "cm")
zoom_width <- grid::unit(config$zoom$width_cm, "cm")
block_names <- character()
label_map <- stats::setNames(shown, shown)

for (db in dbs) {
    block_id <- paste0(".", db, "_block")
    db_tbl <- enr[enr$database == db, , drop = FALSE]
    right_args[[block_id]] <- ComplexHeatmap::anno_block(
        gp = grid::gpar(fill = unname(direct_cols[levels(direct)])),
        which = "row",
        width = grid::unit(3, "mm")
    )
    right_args[[db]] <- ComplexHeatmap::anno_zoom(
        align_to = direct,
        which = "row",
        panel_fun = make_zoom(db, db_tbl),
        size = zoom_size,
        width = zoom_width
    )
    block_names <- c(block_names, block_id)
    label_map[[block_id]] <- ""
    label_map[[db]] <- db
}

anno_ids <- names(right_args)
n_right <- length(right_args)
gap_vec <- rep(1, max(0L, n_right - 1L))
# 0 mm between each direction block and its zoom panel.
if (length(block_names)) {
    for (bn in block_names) {
        i <- match(bn, anno_ids)
        if (!is.na(i) && i < n_right) {
            gap_vec[[i]] <- 0
        }
    }
}

right_args$annotation_name_side <- "top"
right_args$annotation_label <- unname(label_map[anno_ids])
if (n_right > 1L) {
    right_args$gap <- grid::unit(gap_vec, "mm")
}
right_anno <- do.call(ComplexHeatmap::rowAnnotation, right_args)

lgd_list <- list(
    ComplexHeatmap::Legend(
        col_fun = col_fun,
        title = config$labels$fill,
        at = c(-2, -1, 0, 1, 2)
    ),
    ComplexHeatmap::Legend(
        labels = names(direct_cols),
        legend_gp = grid::gpar(fill = unname(direct_cols)),
        title = config$labels$direct
    )
)
if (!is.null(tf_vec) && any(!is.na(tf_vec))) {
    used <- intersect(names(tf_cols), unique(tf_vec[!is.na(tf_vec)]))
    lgd_list <- c(lgd_list, list(
        ComplexHeatmap::Legend(
            labels = used,
            legend_gp = grid::gpar(fill = unname(tf_cols[used])),
            title = config$labels$tf
        )
    ))
}
if (!is.null(sp_vec) && any(!is.na(sp_vec))) {
    used <- intersect(names(sp_cols), unique(sp_vec[!is.na(sp_vec)]))
    lgd_list <- c(lgd_list, list(
        ComplexHeatmap::Legend(
            labels = used,
            legend_gp = grid::gpar(fill = unname(sp_cols[used])),
            title = config$labels$surface
        )
    ))
}
p_at <- pretty(c(0, maxP), n = 4L)
lgd_list <- c(lgd_list, list(
    ComplexHeatmap::Legend(
        at = p_at,
        col_fun = gradient_up,
        title = config$labels$p_up
    ),
    ComplexHeatmap::Legend(
        at = p_at,
        col_fun = gradient_down,
        title = config$labels$p_down
    )
))
lgd <- ComplexHeatmap::packLegend(list = lgd_list)

ht <- ComplexHeatmap::Heatmap(
    mat,
    name = config$labels$fill,
    col = col_fun,
    na_col = "transparent",
    show_heatmap_legend = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    border = TRUE,
    left_annotation = left_anno,
    right_annotation = right_anno,
    column_split = group,
    row_split = direct,
    row_title = NULL,
    use_raster = TRUE
)

# -----------------------------------------------------------------------------
# SAVE
# Nested ggplot grobs in anno_zoom cannot be captured by grid.grabExpr
# (later panels overwrite earlier ones). Draw on the output device.
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
width <- config$size$width
height <- config$size$height
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
        merge_legend = FALSE
    ),
    finally = grDevices::dev.off()
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
