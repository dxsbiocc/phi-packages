#!/usr/bin/env Rscript

# Template-ID: scatter-variate-panel
#
# Purpose:
#   Faceted scatter of supplied component scores (X-variate 1 vs 2).
#   Each panel is one block (for example Human / Microbiome) crossed
#   with one status (stable / progressing). Time is colour and shape
#   together. Supplied max.dist class polygons paint the background
#   of each panel. Axis titles carry that panel's explained-variance
#   percent.
#
# Inputs:
#   1. Scores, one row per sample. Default columns:
#        - label, x, y, time, block, status, expl_x, expl_y
#      Optional label_x, label_y place a sample name by hand.
#      Missing label coordinates fall back to repelled text.
#   2. Background polygons, one row per vertex:
#        - block, status, class, region, ord, x, y
#      `class` is the month predicted inside that polygon. Vertices of
#      one polygon share `region` and follow `ord`. Build both tables
#      with prepare.R from a feature matrix. This script only draws them.
#
# Output:
#   A PDF, PNG, or SVG multi-panel score plot with one shared legend.
#
# Dependencies:
#   ggplot2, readr, ggrepel, patchwork, jsonlite
#
# Example:
#   Rscript prepare.R --example
#   Rscript plot.R example.tsv background.tsv output.pdf
#
# Real data:
#   Rscript prepare.R X.tsv samples.tsv scores.tsv background.tsv
#   Rscript plot.R scores.tsv background.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (column names, time colours and shapes, strip
#   fills, axis-title templates, label toggle). Edit DATA PREPARATION
#   when the table layout differs. Edit PLOT only when the geometry
#   must change. Fit the model in prepare.R, not here.
#
# Scientific assumptions:
#   Scores, explained-variance percents, and the class polygons are
#   supplied. Each polygon is one max.dist class of a 2-component
#   (s)PLS-DA: the class whose predicted dummy score is largest.
#   That score is linear in the two variates, so each class is a
#   convex polygon. This script does not fit that model and does not
#   tessellate points.

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

io <- parse_io_args(n_input = 2L, input_names = c("scores", "background"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        scores = list(
            label = "label",
            x = "x",
            y = "y",
            time = "time",
            block = "block",
            status = "status",
            expl_x = "expl_x",
            expl_y = "expl_y",
            label_x = "label_x",
            label_y = "label_y"
        ),
        background = list(
            block = "block",
            status = "status",
            class = "class",
            region = "region",
            ord = "ord",
            x = "x",
            y = "y"
        )
    ),
    # NULL = first-appearance order in the scores table.
    block_levels = c("Human", "Microbiome"),
    status_levels = c("stable", "progressing"),
    time_levels = c(
        "0 months", "2 months", "4 months", "6 months",
        "8 months", "10 months", "12 months"
    ),
    # Hex copied from colors.json. Unlisted levels fall back to
    # Qualitative.Bold (colours) or default_shapes (shapes).
    time_colours = c(
        "0 months" = "#3969ac",  # Qualitative.Bold
        "2 months" = "#11a579",  # Qualitative.Bold
        "4 months" = "#e73f74",  # Qualitative.Bold
        "6 months" = "#f2b701",  # Qualitative.Bold
        "8 months" = "#000000",  # Brand.adidas
        "10 months" = "#e68310", # Qualitative.Bold
        "12 months" = "#000000"  # Brand.adidas
    ),
    time_shapes = c(
        "0 months" = 16,  # filled circle
        "2 months" = 4,   # x
        "4 months" = 18,  # filled diamond
        "6 months" = 25,  # filled down-triangle
        "8 months" = 7,   # square with x
        "10 months" = 17, # filled triangle
        "12 months" = 3   # plus
    ),
    default_shapes = c(16, 17, 15, 18, 4, 3, 8, 25, 7, 10, 9, 13),
    # Strip fills by status. Qualitative.Pastel1.
    strip_fill = c(
        stable = "#decbe4",
        progressing = "#fbb4ae"
    ),
    axis = list(
        x = "X-variate 1: %s%% expl. var",
        y = "X-variate 2: %s%% expl. var"
    ),
    panel_sep = " ",
    show_labels = TRUE,
    background_alpha = 0.4,
    point_size = 3.3,
    point_stroke = 0.35,
    label_size = 2.55,
    legend_key_size = 3.5
)

load_packages(c("ggplot2", "readr", "ggrepel", "patchwork", "jsonlite"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
resolve_levels <- function(configured, observed) {
    observed <- unique(as.character(observed))
    observed <- observed[!is.na(observed) & nzchar(observed)]
    if (is.null(configured) || !length(configured)) {
        return(observed)
    }
    configured <- as.character(configured)
    c(intersect(configured, observed), setdiff(observed, configured))
}

format_expl <- function(x) {
    x <- as.numeric(x)[[1]]
    if (!is.finite(x)) {
        stop("Explained variance is not finite.", call. = FALSE)
    }
    if (abs(x - round(x)) < 1e-6) {
        sprintf("%d", as.integer(round(x)))
    } else {
        sprintf("%.1f", x)
    }
}

named_or_default <- function(levels, named, fallback) {
    out <- setNames(fallback[seq_along(levels)], levels)
    if (!is.null(named) && length(named)) {
        hit <- intersect(levels, names(named))
        out[hit] <- named[hit]
    }
    out
}

score_cols <- config$columns$scores
bg_cols <- config$columns$background

scores <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
background <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
require_columns(scores, score_cols[c(
    "label", "x", "y", "time", "block", "status", "expl_x", "expl_y"
)])
require_columns(background, bg_cols)

for (nm in c("x", "y", "expl_x", "expl_y")) {
    scores[[score_cols[[nm]]]] <- as.numeric(scores[[score_cols[[nm]]]])
}
background[[bg_cols$x]] <- as.numeric(background[[bg_cols$x]])
background[[bg_cols$y]] <- as.numeric(background[[bg_cols$y]])
background[[bg_cols$ord]] <- as.integer(background[[bg_cols$ord]])

chr_cols <- c("label", "time", "block", "status")
for (nm in chr_cols) {
    scores[[score_cols[[nm]]]] <- as.character(scores[[score_cols[[nm]]]])
}
for (nm in c("block", "status", "class", "region")) {
    background[[bg_cols[[nm]]]] <- as.character(background[[bg_cols[[nm]]]])
}

score_ok <- is.finite(scores[[score_cols$x]]) &
    is.finite(scores[[score_cols$y]]) &
    is.finite(scores[[score_cols$expl_x]]) &
    is.finite(scores[[score_cols$expl_y]]) &
    !is.na(scores[[score_cols$time]]) & nzchar(scores[[score_cols$time]]) &
    !is.na(scores[[score_cols$block]]) & nzchar(scores[[score_cols$block]]) &
    !is.na(scores[[score_cols$status]]) & nzchar(scores[[score_cols$status]])
if (any(!score_ok)) {
    message("Dropped ", sum(!score_ok), " score row(s) with missing coordinates, variance, time, block, or status.")
}
scores <- scores[score_ok, , drop = FALSE]
if (!nrow(scores)) {
    stop("No score rows left to plot.", call. = FALSE)
}
if (any(is.na(scores[[score_cols$label]]))) {
    scores[[score_cols$label]][is.na(scores[[score_cols$label]])] <- ""
}

bg_ok <- is.finite(background[[bg_cols$x]]) &
    is.finite(background[[bg_cols$y]]) &
    is.finite(background[[bg_cols$ord]]) &
    !is.na(background[[bg_cols$class]]) & nzchar(background[[bg_cols$class]]) &
    !is.na(background[[bg_cols$region]]) & nzchar(background[[bg_cols$region]]) &
    !is.na(background[[bg_cols$block]]) & nzchar(background[[bg_cols$block]]) &
    !is.na(background[[bg_cols$status]]) & nzchar(background[[bg_cols$status]])
if (any(!bg_ok)) {
    message("Dropped ", sum(!bg_ok), " background row(s) with non-finite coordinates or empty class, region, block, or status.")
}
background <- background[bg_ok, , drop = FALSE]
if (!nrow(background)) {
    stop("No background vertices left to plot.", call. = FALSE)
}

expl_x <- scores[[score_cols$expl_x]]
expl_y <- scores[[score_cols$expl_y]]
if (any(expl_x < 0 | expl_x > 100 | expl_y < 0 | expl_y > 100)) {
    stop(
        "expl_x and expl_y must be percents in [0, 100], not proportions.",
        call. = FALSE
    )
}

block_levels <- resolve_levels(config$block_levels, scores[[score_cols$block]])
status_levels <- resolve_levels(config$status_levels, scores[[score_cols$status]])
time_levels <- resolve_levels(config$time_levels, scores[[score_cols$time]])
if (length(block_levels) < 1L || length(status_levels) < 1L) {
    stop("Need at least one block and one status.", call. = FALSE)
}
if (length(time_levels) < 2L) {
    stop("Need at least two time levels so colour and shape encode time.", call. = FALSE)
}

extra_class <- setdiff(unique(background[[bg_cols$class]]), time_levels)
if (length(extra_class)) {
    stop(
        paste(
            "Background class is not a time level:",
            paste(extra_class, collapse = ", "),
            "\nUse the same labels as the scores time column."
        ),
        call. = FALSE
    )
}

scores[[score_cols$block]] <- factor(scores[[score_cols$block]], levels = block_levels)
scores[[score_cols$status]] <- factor(scores[[score_cols$status]], levels = status_levels)
scores[[score_cols$time]] <- factor(scores[[score_cols$time]], levels = time_levels)
background[[bg_cols$block]] <- factor(background[[bg_cols$block]], levels = block_levels)
background[[bg_cols$status]] <- factor(background[[bg_cols$status]], levels = status_levels)
background[[bg_cols$class]] <- factor(background[[bg_cols$class]], levels = time_levels)

expected <- expand.grid(
    block = block_levels,
    status = status_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
present <- unique(data.frame(
    block = as.character(scores[[score_cols$block]]),
    status = as.character(scores[[score_cols$status]]),
    stringsAsFactors = FALSE
))
missing_panels <- expected[
    !paste(expected$block, expected$status) %in% paste(present$block, present$status),
    ,
    drop = FALSE
]
if (nrow(missing_panels)) {
    stop(
        paste(
            "Every block x status pair needs scores. Missing:",
            paste(
                paste(missing_panels$block, missing_panels$status),
                collapse = "; "
            )
        ),
        call. = FALSE
    )
}

bg_present <- unique(data.frame(
    block = as.character(background[[bg_cols$block]]),
    status = as.character(background[[bg_cols$status]]),
    stringsAsFactors = FALSE
))
missing_bg <- expected[
    !paste(expected$block, expected$status) %in% paste(bg_present$block, bg_present$status),
    ,
    drop = FALSE
]
if (nrow(missing_bg)) {
    stop(
        paste(
            "Every block x status pair needs class polygons. Missing:",
            paste(paste(missing_bg$block, missing_bg$status), collapse = "; ")
        ),
        call. = FALSE
    )
}

variance_ok <- function(block, status) {
    keep <- scores[[score_cols$block]] == block &
        scores[[score_cols$status]] == status
    vx <- unique(scores[[score_cols$expl_x]][keep])
    vy <- unique(scores[[score_cols$expl_y]][keep])
    if (length(vx) != 1L || length(vy) != 1L) {
        stop(
            paste0(
                "expl_x and expl_y must be constant within ",
                block, config$panel_sep, status, "."
            ),
            call. = FALSE
        )
    }
    c(vx, vy)
}

colour_fallback <- palette_colors("Qualitative.Bold", n = length(time_levels))
time_colours <- named_or_default(time_levels, config$time_colours, colour_fallback)
# Polygons are drawn opaque. Mix the class colour toward white here so
# shared edges do not leave the seams an alpha tile would.
lighten <- function(hex, alpha) {
    rgb <- grDevices::col2rgb(hex) / 255
    mixed <- rgb * alpha + (1 - alpha)
    grDevices::rgb(mixed[1, ], mixed[2, ], mixed[3, ])
}
region_fill <- setNames(
    lighten(unname(time_colours), as.numeric(config$background_alpha)),
    names(time_colours)
)

shape_pool <- as.integer(config$default_shapes)
if (!length(shape_pool)) {
    stop("default_shapes must contain at least one ggplot shape.", call. = FALSE)
}
shape_fallback <- shape_pool[((seq_along(time_levels) - 1L) %% length(shape_pool)) + 1L]
time_shapes <- as.integer(named_or_default(
    time_levels, config$time_shapes, shape_fallback
))
if (any(!is.finite(time_shapes) | time_shapes < 0 | time_shapes > 25)) {
    stop("time_shapes must be ggplot shape numbers from 0 to 25.", call. = FALSE)
}

strip_fallback <- palette_colors("Qualitative.Pastel1", n = length(status_levels))
strip_fill <- named_or_default(status_levels, config$strip_fill, strip_fallback)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
make_panel <- function(block, status) {
    pts <- scores[
        scores[[score_cols$block]] == block &
            scores[[score_cols$status]] == status,
        ,
        drop = FALSE
    ]
    bg <- background[
        background[[bg_cols$block]] == block &
            background[[bg_cols$status]] == status,
        ,
        drop = FALSE
    ]
    bg <- bg[order(bg[[bg_cols$region]], bg[[bg_cols$ord]]), , drop = FALSE]
    region_n <- tapply(bg[[bg_cols$ord]], bg[[bg_cols$region]], length)
    if (any(region_n < 3L)) {
        stop(paste("A class polygon in", block, status, "has fewer than 3 vertices."), call. = FALSE)
    }
    region_dup <- tapply(
        bg[[bg_cols$ord]],
        bg[[bg_cols$region]],
        function(v) anyDuplicated(v) > 0
    )
    if (any(unlist(region_dup))) {
        stop(paste("Polygon vertex order repeats in", block, status), call. = FALSE)
    }

    expl <- variance_ok(block, status)
    x_title <- sprintf(config$axis$x, format_expl(expl[[1]]))
    y_title <- sprintf(config$axis$y, format_expl(expl[[2]]))
    panel_title <- paste(block, status, sep = config$panel_sep)
    pts$.panel <- factor(panel_title)
    bg$.panel <- factor(panel_title)
    bg$.fill <- unname(region_fill[as.character(bg[[bg_cols$class]])])
    pts$.fill <- unname(time_colours[as.character(pts[[score_cols$time]])])
    xlim <- range(c(bg[[bg_cols$x]], pts[[score_cols$x]]))
    ylim <- range(c(bg[[bg_cols$y]], pts[[score_cols$y]]))

    p <- ggplot() +
        geom_polygon(
            data = bg,
            aes(
                x = .data[[bg_cols$x]],
                y = .data[[bg_cols$y]],
                group = .data[[bg_cols$region]],
                fill = .data$.fill,
                colour = after_scale(fill)
            ),
            linewidth = 0.4,
            show.legend = FALSE
        ) +
        geom_point(
            data = pts,
            aes(
                x = .data[[score_cols$x]],
                y = .data[[score_cols$y]],
                colour = .data[[score_cols$time]],
                fill = .data$.fill,
                shape = .data[[score_cols$time]]
            ),
            size = as.numeric(config$point_size),
            stroke = as.numeric(config$point_stroke)
        ) +
        facet_wrap(vars(.panel)) +
        scale_x_continuous(
            name = x_title,
            limits = xlim,
            expand = c(0, 0)
        ) +
        scale_y_continuous(
            name = y_title,
            limits = ylim,
            expand = c(0, 0)
        ) +
        scale_colour_manual(
            values = time_colours,
            breaks = time_levels,
            drop = FALSE,
            name = NULL
        ) +
        scale_fill_identity(guide = "none") +
        scale_shape_manual(
            values = time_shapes,
            breaks = time_levels,
            drop = FALSE,
            name = NULL
        ) +
        guides(
            colour = guide_legend(
                order = 1,
                override.aes = list(
                    size = as.numeric(config$legend_key_size),
                    alpha = 1,
                    fill = unname(time_colours),
                    colour = unname(time_colours)
                )
            ),
            shape = guide_legend(order = 1)
        ) +
        theme_classic(base_size = 10) +
        theme(
            aspect.ratio = 1,
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            strip.background = element_rect(
                fill = unname(strip_fill[[status]]),
                colour = NA
            ),
            strip.text = element_text(
                face = "bold",
                size = 12.5,
                colour = "#222222",
                margin = margin(5, 4, 5, 4)
            ),
            axis.title = element_text(size = 9, colour = "#333333"),
            axis.text = element_text(size = 8, colour = "#4a4a4a"),
            axis.line = element_line(colour = "#4a4a4a", linewidth = 0.35),
            axis.ticks = element_line(colour = "#4a4a4a", linewidth = 0.3),
            legend.position = "left",
            legend.key = element_blank(),
            legend.background = element_blank(),
            legend.text = element_text(size = 10.5, colour = "#222222"),
            legend.key.height = unit(15, "pt"),
            legend.key.width = unit(16, "pt"),
            legend.spacing.y = unit(1.5, "pt"),
            plot.margin = margin(3, 8, 3, 6)
        )

    if (isTRUE(config$show_labels)) {
        lab <- pts[nzchar(pts[[score_cols$label]]), , drop = FALSE]
        lx_col <- score_cols$label_x
        ly_col <- score_cols$label_y
        fixed <- rep(FALSE, nrow(lab))
        if (nrow(lab) && !is.null(lx_col) && !is.null(ly_col) &&
            all(c(lx_col, ly_col) %in% names(lab))) {
            lab[[lx_col]] <- as.numeric(lab[[lx_col]])
            lab[[ly_col]] <- as.numeric(lab[[ly_col]])
            fixed <- is.finite(lab[[lx_col]]) & is.finite(lab[[ly_col]])
        }
        fixed_lab <- lab[fixed, , drop = FALSE]
        free_lab <- lab[!fixed, , drop = FALSE]
        if (nrow(fixed_lab)) {
            p <- p +
                geom_text(
                    data = fixed_lab,
                    aes(
                        x = .data[[lx_col]],
                        y = .data[[ly_col]],
                        label = .data[[score_cols$label]]
                    ),
                    inherit.aes = FALSE,
                    size = as.numeric(config$label_size),
                    colour = "#222222",
                    hjust = 0.5,
                    vjust = 0.5
                )
        }
        if (nrow(free_lab)) {
            p <- p +
                ggrepel::geom_text_repel(
                    data = free_lab,
                    aes(
                        x = .data[[score_cols$x]],
                        y = .data[[score_cols$y]],
                        label = .data[[score_cols$label]]
                    ),
                    size = as.numeric(config$label_size),
                    colour = "#222222",
                    segment.colour = "#888888",
                    segment.size = 0.25,
                    min.segment.length = 0.2,
                    box.padding = 0.28,
                    point.padding = 0.18,
                    force = 2.4,
                    force_pull = 0.35,
                    max.overlaps = Inf,
                    max.time = 2,
                    xlim = xlim,
                    ylim = ylim,
                    seed = 1,
                    show.legend = FALSE
                )
        }
    }
    p
}

panels <- list()
for (block in block_levels) {
    for (status in status_levels) {
        panels[[length(panels) + 1L]] <- make_panel(block, status)
    }
}

# Collected legend sits in its own left column so it is shared across
# the block x status grid. guide_area() avoids patchwork's `&` theme
# helper, which this ggplot2 / patchwork pair does not dispatch.
if (length(panels) > 25L) {
    stop(
        "This layout supports at most 25 panels (block x status).",
        call. = FALSE
    )
}
if (length(block_levels) == 1L && length(status_levels) == 1L) {
    p <- panels[[1]]
} else {
    design <- paste(
        vapply(seq_along(block_levels), function(i) {
            start <- (i - 1L) * length(status_levels) + 1L
            cells <- vapply(
                seq_len(length(status_levels)),
                function(j) LETTERS[start + j - 1L],
                character(1)
            )
            paste0("#", paste(cells, collapse = ""))
        }, character(1)),
        collapse = "\n"
    )
    # The legend column is marked '#' only as a placeholder; the real
    # area is the guide plot appended last, lettered after the panels.
    legend_letter <- LETTERS[length(panels) + 1L]
    design <- gsub("#", legend_letter, design, fixed = TRUE)
    p <- patchwork::wrap_plots(c(panels, list(patchwork::guide_area()))) +
        patchwork::plot_layout(
            design = design,
            guides = "collect",
            widths = c(0.62, rep(1, length(status_levels)))
        )
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
n_status <- length(status_levels)
n_block <- length(block_levels)
save_ggplot(
    p,
    io$output,
    width = max(7.6, 1.85 + 3.65 * n_status),
    height = max(6.2, 3.85 * n_block)
)
