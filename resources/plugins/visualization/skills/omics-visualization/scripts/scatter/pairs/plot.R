#!/usr/bin/env Rscript

# Template-ID: scatter-pairs
#
# Purpose:
#   Draw a base-graphics pairs grid of several numeric variables:
#   lower triangle scatter plus lm line, diagonal variable names,
#   upper triangle Pearson r as coloured tiles with printed
#   coefficients. Matches the classic pairs custom-panel look.
#
# Inputs:
#   One row per observation. Default example:
#     - sample: optional id (dropped from the grid)
#     - remaining numeric columns: variables (genes)
#
# Output:
#   A PDF, PNG, or SVG pairs plot.
#
# Dependencies:
#   readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Extra variables are extra numeric columns, not extra ids.
#   Correlation method and gap stay in CONFIG. One pair is
#   scatter-correlation. Coefficients only, no point clouds, is
#   heatmap-signif. Category × category bubbles are scatter-matrix.
#   Do not run DE, enrichment, or clustering here.
#
# Scientific assumptions:
#   Each row is one paired observation of the plotted variables.
#   r is Pearson (or CONFIG method) of the displayed points with
#   pairwise complete observations. The lm line is OLS, not a
#   causal claim. More than about 12 variables is usually unreadable
#   here; use heatmap-signif.

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
    variables = NULL,
    method = "pearson",
    gap = 0.5,
    palettes = list(
        fill = "Diverging.Temps",
        point = "Brand.Algolia",
        line = "Brand.Algolia"
    ),
    labels = list(
        title = ""
    ),
    size = list(
        width = NULL,
        height = NULL
    )
)

load_packages(c("readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
id_col <- config$columns$id
if (!is.null(id_col) && nzchar(id_col) && id_col %in% names(df)) {
    df[[id_col]] <- NULL
}

if (is.null(config$variables)) {
    vars <- names(df)
} else {
    vars <- as.character(config$variables)
    missing <- setdiff(vars, names(df))
    if (length(missing)) {
        stop(
            paste("Missing variable column(s):", paste(missing, collapse = ", ")),
            call. = FALSE
        )
    }
}

for (nm in vars) {
    df[[nm]] <- as.numeric(df[[nm]])
}
keep <- vapply(vars, function(nm) {
    any(is.finite(df[[nm]]))
}, logical(1))
vars <- vars[keep]
if (length(vars) < 2L) {
    stop("Need at least two numeric variable columns.", call. = FALSE)
}
mat <- as.data.frame(df[, vars, drop = FALSE], optional = TRUE)

n_var <- length(vars)
if (n_var > 12L) {
    message(
        n_var, " variables: pairs cells get small. heatmap-signif is",
        " usually clearer for a coefficient-only matrix."
    )
}

brand <- palette_colors(config$palettes$point)
point_col <- grDevices::adjustcolor(brand[[1]], alpha.f = 0.63)
line_col <- palette_colors(config$palettes$line)[[5]]
fill_stops <- palette_colors(config$palettes$fill)
method <- config$method
gap <- as.numeric(config$gap)
if (!is.finite(gap) || gap < 0) {
    gap <- 0.5
}
n_fill <- 14L
fill_pal <- grDevices::colorRampPalette(fill_stops)(n_fill)
outline <- "black"

panel_lower <- function(x, y, corr = NULL, ...) {
    if (!is.null(corr)) {
        return(invisible())
    }
    graphics::plot.xy(
        grDevices::xy.coords(x, y),
        type = "p",
        pch = 19,
        cex = 0.8,
        col = point_col,
        ...
    )
    fit <- stats::lm(y ~ x)
    graphics::abline(fit, lwd = 2, col = line_col)
}

panel_diag <- function(x = 0.5, y = 0.5, txt, cex, font) {
    graphics::text(x, y, txt, cex = cex, font = font)
    graphics::box(col = outline, lwd = 2)
}

panel_upper <- function(x, y, corr = NULL, ...) {
    r <- suppressWarnings(
        stats::cor(x, y, use = "pairwise.complete.obs", method = method)
    )
    if (!is.finite(r)) {
        r <- 0
    }
    col_ind <- as.numeric(cut(
        r,
        breaks = seq(from = -1, to = 1, length.out = n_fill + 1L),
        include.lowest = TRUE
    ))
    graphics::par(new = TRUE)
    graphics::plot(
        0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1),
        axes = FALSE, asp = 1
    )
    usr <- graphics::par("usr")
    graphics::rect(
        usr[[1]], usr[[3]], usr[[2]], usr[[4]],
        col = fill_pal[[col_ind]],
        border = NA
    )
    txt_col <- if (.hex_luminance(fill_pal[[col_ind]]) > 0.55) {
        "black"
    } else {
        "white"
    }
    graphics::text(
        0, 0,
        labels = sprintf("%.2f", r),
        cex = 2.2,
        col = txt_col,
        font = 2
    )
    graphics::box(col = outline, lwd = 1)
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
width <- config$size$width
height <- config$size$height
if (is.null(width) || !is.finite(as.numeric(width))) {
    width <- max(8, 0.72 * n_var + 2.2)
}
if (is.null(height) || !is.finite(as.numeric(height))) {
    height <- width
}
width <- as.numeric(width)
height <- as.numeric(height)
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
op <- graphics::par(no.readonly = TRUE)
tryCatch({
    graphics::par(bg = NA, mar = c(2.2, 2.2, 2.2, 2.2))
    if (nzchar(config$labels$title)) {
        graphics::par(oma = c(0, 0, 2, 0))
    }
    graphics::pairs(
        mat,
        gap = gap,
        text.panel = panel_diag,
        lower.panel = panel_lower,
        upper.panel = panel_upper,
        cex.labels = 1.35,
        font.labels = 2
    )
    if (nzchar(config$labels$title)) {
        graphics::mtext(
            config$labels$title,
            side = 3,
            outer = TRUE,
            line = 0.4,
            cex = 1.1,
            font = 2
        )
    }
}, finally = {
    graphics::par(op)
    grDevices::dev.off()
})
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
