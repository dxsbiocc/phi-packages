#!/usr/bin/env Rscript

# Template-ID: scatter-bubble
#
# Purpose:
#   Draw a bubble scatter of two numeric axes with point size from a
#   third magnitude, colour by a discrete group, optional log10 x-axis,
#   optional overall linear fit, selective id labels, optional nested
#   size key, and optional cluster callouts (preview.png;
#   preview_methane.png).
#
# Inputs:
#   One row per observation. Default example:
#     - gni, coverage, fractions, income_group, country, show_label
#   Methane variant example_methane.tsv:
#     - diff_share, invent_share, emission_share, country_group, code,
#       show_label
#   Optional callouts TSV (callouts_path):
#     - label, x, y, xend, yend
#
# Output:
#   A PDF, PNG, or SVG bubble chart.
#
# Dependencies:
#   ggplot2, readr, scales, ggrepel, grid
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, colours, size guide, fit, labels, limits,
#   and callouts. Extra groups or labelled ids are extra rows, not a
#   new id. A dual-panel threshold bubble stays scatter-bubble-panel.
#   A category × category bubble grid stays scatter-matrix.
#
# Scientific assumptions:
#   x, y, size, and group are supplied. This script does not compute
#   GNI, coverage, invention/emission shares, or share differences.
#   The optional lm line is a display fit of y ~ log10(x) when
#   x_log10 is TRUE (else y ~ x). Bubble area encodes the size column.

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
        x = "gni",
        y = "coverage",
        size = "fractions",
        group = "income_group",
        id = "country",
        show_label = "show_label"
    ),
    group_order = c(
        "High-income countries",
        "Upper-middle-income countries",
        "Lower-middle-income countries",
        "Low-income countries"
    ),
    group_colours = c(
        `High-income countries` = "#00A8B5",
        `Upper-middle-income countries` = "#E07A3D",
        `Lower-middle-income countries` = "#5BA33B",
        `Low-income countries` = "#A51C4A"
    ),
    palettes = list(group = "Qualitative.Safe"),
    x_log10 = TRUE,
    x_breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000),
    xlim = c(180, 120000),
    ylim = c(0, 1.55),
    # Size legend reference values (in data units of the size column).
    size_breaks = c(1e5, 2e6, 4e6, 6e6, 8e6, 1e7),
    size_labels = c(
        "<= 100,000", "2,000,000", "4,000,000",
        "6,000,000", "8,000,000", ">= 10,000,000"
    ),
    # "legend" = ggplot size guide (with colour);
    # "nested" = concentric key on the right;
    # "right" = stacked discrete size key on the right.
    size_guide = "legend",
    size_max = 14,
    nested_size_x = NULL,
    nested_size_y = NULL,
    # For size_guide = "right": anchor of the stacked key (data coords).
    right_size_x = NULL,
    right_size_y = NULL,
    point_stroke = 0.85,
    point_alpha = 0.95,
    # Hollow open circles (shape 1). Set FALSE for filled 16.
    hollow = TRUE,
    show_fit = TRUE,
    fit_colour = "grey35",
    fit_linetype = "dashed",
    fit_linewidth = 0.55,
    # Label every id when TRUE (ignores show_label).
    label_all = FALSE,
    label_size = 3.1,
    label_nudge_y = 0.02,
    # Optional callouts TSV beside the script, or absolute path; NULL = none.
    callouts_path = NULL,
    seed = 1L,
    labels = list(
        title = NULL,
        x = "Gross national income, 2013 (US$)",
        y = "Coverage",
        size = "Fractions per year",
        group = "Income group"
    ),
    legend_position = c(0.02, 0.98),
    size = list(width = 8.2, height = 6.4)
)

load_packages(c("ggplot2", "readr", "scales", "ggrepel", "grid"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
req <- config$columns[c("x", "y", "size", "group")]
require_columns(df, req)

x_col <- config$columns$x
y_col <- config$columns$y
sz_col <- config$columns$size
g_col <- config$columns$group
id_col <- config$columns$id
lab_flag <- config$columns$show_label

df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df[[sz_col]] <- as.numeric(df[[sz_col]])
df[[g_col]] <- as.character(df[[g_col]])

ok <- is.finite(df[[x_col]]) &
    is.finite(df[[y_col]]) & is.finite(df[[sz_col]]) & df[[sz_col]] >= 0 &
    !is.na(df[[g_col]]) & nzchar(df[[g_col]])
if (isTRUE(config$x_log10)) {
    ok <- ok & df[[x_col]] > 0
}
df <- df[ok, , drop = FALSE]
if (!nrow(df)) stop("No finite bubble rows to plot.", call. = FALSE)

g_lv <- if (!is.null(config$group_order)) {
    c(
        intersect(as.character(config$group_order), unique(df[[g_col]])),
        setdiff(unique(df[[g_col]]), as.character(config$group_order))
    )
} else {
    unique(df[[g_col]])
}
df[[g_col]] <- factor(df[[g_col]], levels = g_lv)

g_cols <- config$group_colours
g_cols <- g_cols[names(g_cols) %in% g_lv]
miss <- setdiff(g_lv, names(g_cols))
if (length(miss)) {
    extra <- expand_palette(palette_colors(config$palettes$group), length(miss))
    names(extra) <- miss
    g_cols <- c(g_cols, extra)
}
g_cols <- g_cols[g_lv]

df$`.label` <- ""
if (!is.null(id_col) && id_col %in% names(df)) {
    df[[id_col]] <- as.character(df[[id_col]])
    if (isTRUE(config$label_all)) {
        df$`.label` <- df[[id_col]]
    } else if (!is.null(lab_flag) && lab_flag %in% names(df)) {
        show <- as.integer(df[[lab_flag]])
        show[is.na(show)] <- 0L
        df$`.label`[show > 0L] <- df[[id_col]][show > 0L]
    }
}
lab <- df[nzchar(df$`.label`), , drop = FALSE]

sz_breaks <- as.numeric(config$size_breaks)
sz_labs <- config$size_labels
if (is.null(sz_labs) || length(sz_labs) != length(sz_breaks)) {
    sz_labs <- scales::comma(sz_breaks)
}

# Optional cluster callouts.
callouts <- NULL
cpath <- config$callouts_path
if (!is.null(cpath) && nzchar(cpath)) {
    if (!grepl("^(/|[A-Za-z]:)", cpath)) {
        script_dir <- dirname(normalizePath(sub(
            "^--file=", "",
            grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]]
        )))
        cpath <- file.path(script_dir, cpath)
    }
    if (file.exists(cpath)) {
        callouts <- as.data.frame(read_table_auto(cpath), stringsAsFactors = FALSE)
        need_c <- c("label", "x", "y", "xend", "yend")
        if (!all(need_c %in% names(callouts))) {
            stop(
                "callouts TSV needs columns: label, x, y, xend, yend",
                call. = FALSE
            )
        }
        callouts$x <- as.numeric(callouts$x)
        callouts$y <- as.numeric(callouts$y)
        callouts$xend <- as.numeric(callouts$xend)
        callouts$yend <- as.numeric(callouts$yend)
        callouts$label <- gsub("\\\\n", "\n", as.character(callouts$label))
        callouts <- callouts[
            is.finite(callouts$x) & is.finite(callouts$y) &
                is.finite(callouts$xend) & is.finite(callouts$yend),
            ,
            drop = FALSE
        ]
    }
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
pt_shape <- if (isTRUE(config$hollow)) 1L else 16L
size_guide_mode <- as.character(config$size_guide)
use_nested <- identical(size_guide_mode, "nested")
use_right_size <- identical(size_guide_mode, "right")
hide_size_guide <- use_nested || use_right_size

p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]]))

if (isTRUE(config$show_fit)) {
    if (isTRUE(config$x_log10)) {
        fit_df <- data.frame(
            x = df[[x_col]],
            y = df[[y_col]],
            stringsAsFactors = FALSE
        )
        fit_df <- fit_df[is.finite(fit_df$x) & fit_df$x > 0 & is.finite(fit_df$y), ]
        mod <- stats::lm(y ~ log10(x), data = fit_df)
        xr <- range(fit_df$x, finite = TRUE)
        xs <- 10^seq(log10(xr[[1]]), log10(xr[[2]]), length.out = 80)
        fit_line <- data.frame(
            x = xs,
            y = stats::predict(mod, newdata = data.frame(x = xs)),
            stringsAsFactors = FALSE
        )
        p <- p +
            geom_line(
                data = fit_line,
                aes(x = x, y = y),
                inherit.aes = FALSE,
                colour = config$fit_colour,
                linetype = config$fit_linetype,
                linewidth = config$fit_linewidth
            )
    } else {
        p <- p +
            geom_smooth(
                method = "lm",
                formula = y ~ x,
                se = FALSE,
                colour = config$fit_colour,
                linetype = config$fit_linetype,
                linewidth = config$fit_linewidth
            )
    }
}

p <- p +
    geom_point(
        aes(size = .data[[sz_col]], colour = .data[[g_col]]),
        shape = pt_shape,
        stroke = config$point_stroke,
        alpha = config$point_alpha
    )

x_scale <- if (isTRUE(config$x_log10)) {
    scale_x_log10(
        breaks = as.numeric(config$x_breaks),
        labels = scales::label_comma(accuracy = 1)
    )
} else if (!is.null(config$x_breaks)) {
    scale_x_continuous(breaks = as.numeric(config$x_breaks))
} else {
    scale_x_continuous()
}

size_guide_obj <- if (hide_size_guide) {
    "none"
} else if (isTRUE(config$hollow)) {
    guide_legend(
        order = 1,
        override.aes = list(
            colour = "#8AA0B0",
            shape = 1L,
            alpha = 1,
            stroke = 0.7
        )
    )
} else {
    # Filled disks matching solid bubbles (not hollow rings).
    guide_legend(
        order = 1,
        override.aes = list(
            colour = "#6B7C8A",
            shape = 16L,
            alpha = 0.85,
            stroke = 0
        )
    )
}

p <- p +
    x_scale +
    scale_colour_manual(
        name = if (is.null(config$labels$group)) NULL else config$labels$group,
        values = g_cols,
        guide = guide_legend(
            order = 2,
            title = if (is.null(config$labels$group) || !nzchar(config$labels$group)) {
                NULL
            } else {
                config$labels$group
            },
            override.aes = list(
                shape = 16, size = 3.4, alpha = 1, stroke = 0
            )
        )
    ) +
    scale_size_area(
        name = config$labels$size,
        max_size = as.numeric(config$size_max),
        breaks = sz_breaks,
        labels = sz_labs,
        guide = size_guide_obj
    ) +
    coord_cartesian(
        xlim = as.numeric(config$xlim),
        ylim = as.numeric(config$ylim),
        clip = "off"
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_classic(base_size = 12) +
    theme(
        panel.grid = element_blank(),
        axis.line = element_line(colour = "grey20", linewidth = 0.45),
        axis.ticks = element_line(colour = "grey20", linewidth = 0.35),
        legend.position = "inside",
        legend.position.inside = as.numeric(config$legend_position),
        legend.justification = c(0, 1),
        legend.box = "vertical",
        legend.box.just = "left",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 8.5),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.spacing.y = unit(0.15, "lines"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 18, 8, 10)
    )

# Stacked discrete size key on the right (separate circles, not nested).
if (use_right_size) {
    xr <- as.numeric(config$xlim)
    yr <- as.numeric(config$ylim)
    rx <- if (!is.null(config$right_size_x) && is.finite(config$right_size_x)) {
        as.numeric(config$right_size_x)
    } else {
        xr[[1]] + 0.92 * diff(xr)
    }
    ry0 <- if (!is.null(config$right_size_y) && is.finite(config$right_size_y)) {
        as.numeric(config$right_size_y)
    } else {
        yr[[1]] + 0.62 * diff(yr)
    }
    n_sz <- length(sz_breaks)
    # Top = largest; space evenly downward.
    gap <- 0.09 * diff(yr)
    right_df <- data.frame(
        x = rx,
        y = ry0 - (seq_len(n_sz) - 1L) * gap,
        size = rev(sz_breaks),
        lab = rev(sz_labs),
        stringsAsFactors = FALSE
    )
    right_df$lab_x <- rx + 0.035 * diff(xr)
    if (isTRUE(config$hollow)) {
        p <- p +
            geom_point(
                data = right_df,
                aes(x = x, y = y, size = size),
                inherit.aes = FALSE,
                shape = 1L,
                colour = "grey25",
                stroke = 0.55,
                show.legend = FALSE
            )
    } else {
        p <- p +
            geom_point(
                data = right_df,
                aes(x = x, y = y, size = size),
                inherit.aes = FALSE,
                shape = 16L,
                colour = "grey45",
                alpha = 0.85,
                show.legend = FALSE
            )
    }
    p <- p +
        geom_text(
            data = right_df,
            aes(x = lab_x, y = y, label = lab),
            inherit.aes = FALSE,
            hjust = 0,
            size = 2.7,
            colour = "grey20",
            show.legend = FALSE
        )
}

# Nested concentric size key (right side).
if (use_nested) {
    xr <- as.numeric(config$xlim)
    yr <- as.numeric(config$ylim)
    nx <- if (!is.null(config$nested_size_x) && is.finite(config$nested_size_x)) {
        as.numeric(config$nested_size_x)
    } else {
        xr[[1]] + 0.90 * diff(xr)
    }
    ny <- if (!is.null(config$nested_size_y) && is.finite(config$nested_size_y)) {
        as.numeric(config$nested_size_y)
    } else {
        yr[[1]] + 0.52 * diff(yr)
    }
    nest <- data.frame(
        x = nx,
        y = ny,
        size = sz_breaks,
        lab = sz_labs,
        stringsAsFactors = FALSE
    )
    nest <- nest[order(-nest$size), , drop = FALSE]
    n_nest <- nrow(nest)
    nest$lab_x <- nx + 0.07 * diff(xr)
    nest$lab_y <- ny + seq(
        from = 0.10 * diff(yr),
        to = -0.02 * diff(yr),
        length.out = n_nest
    )
    # Match data point style: filled disks vs hollow rings.
    if (isTRUE(config$hollow)) {
        p <- p +
            geom_point(
                data = nest,
                aes(x = x, y = y, size = size),
                inherit.aes = FALSE,
                shape = 1L,
                colour = "grey25",
                stroke = 0.55,
                show.legend = FALSE
            )
    } else {
        p <- p +
            geom_point(
                data = nest,
                aes(x = x, y = y, size = size),
                inherit.aes = FALSE,
                shape = 21L,
                fill = "grey75",
                colour = "grey25",
                stroke = 0.45,
                alpha = 0.9,
                show.legend = FALSE
            )
    }
    p <- p +
        geom_text(
            data = nest,
            aes(x = lab_x, y = lab_y, label = lab),
            inherit.aes = FALSE,
            hjust = 0,
            size = 2.7,
            colour = "grey20",
            show.legend = FALSE
        )
}

# Country labels last: one ggrepel layer only (segments + text stay aligned).
if (nrow(lab)) {
    set.seed(as.integer(config$seed))
    p <- p +
        ggrepel::geom_text_repel(
            data = lab,
            aes(
                x = .data[[x_col]],
                y = .data[[y_col]],
                label = .label,
                colour = .data[[g_col]],
                point.size = .data[[sz_col]]
            ),
            inherit.aes = FALSE,
            size = config$label_size,
            fontface = "bold",
            show.legend = FALSE,
            max.overlaps = Inf,
            force = 1.5,
            force_pull = 0.5,
            box.padding = 0.35,
            point.padding = 0.4,
            min.segment.length = 0.15,
            segment.colour = "grey55",
            segment.size = 0.3,
            seed = as.integer(config$seed),
            xlim = c(-Inf, Inf),
            ylim = c(-Inf, Inf)
        )
}
# Optional cluster callouts (CONFIG callouts_path). Prefer leaving NULL and
# labelling only large points via show_label; callouts are a second style.
if (!is.null(callouts) && nrow(callouts)) {
    p <- p +
        ggrepel::geom_text_repel(
            data = callouts,
            aes(x = xend, y = yend, label = label),
            inherit.aes = FALSE,
            colour = "grey25",
            size = 2.5,
            fontface = "plain",
            show.legend = FALSE,
            max.overlaps = Inf,
            force = 1.2,
            box.padding = 0.3,
            point.padding = 0.2,
            min.segment.length = 0.1,
            segment.colour = "grey50",
            segment.size = 0.3,
            seed = as.integer(config$seed) + 1L,
            xlim = c(-Inf, Inf),
            ylim = c(-Inf, Inf)
        )
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = config$size$width,
    height = config$size$height
)
