#!/usr/bin/env Rscript

# Template-ID: bar-swimmer
#
# Purpose:
#   Draw a clinical swimmer plot: one lane per patient, treatment
#   intervals as rectangles on a time axis, with response points,
#   stop marks, and optional left-hand annotation tiles.
#
# Inputs:
#   One row per treatment interval. Default example:
#     - patient: lane id
#     - treatment: fill of the rectangle (CDK, Other, PARPi, ...)
#     - start, end: time in months (already aligned)
#     - response, response_time: best response and its scan time
#     - stop: Toxicity / Disease Progression / Ongoing
#     - line, loh, cdk, endocrine, platinum: optional annotation columns
#
# Output:
#   A PDF, PNG, or SVG swimmer plot.
#
# Dependencies:
#   ggplot2, ggprism, ggnewscale, patchwork, readr, grid
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, meta, sort, palettes).
#   Extra treatments are extra fill levels. Extra left-hand tiles are
#   extra names in config$meta, not extra ids. Do not add an id per
#   Nature paper layout. Do not compute best response, progression,
#   or Kaplan-Meier curves. A two-point before/after display is
#   scatter-dumbbell. A population survival curve is line-survival.
#
# Scientific assumptions:
#   Start, end, response, and stop are supplied. This script does not
#   convert calendar dates, impute scans, or test treatments. Patient
#   order is a display sort. Arrow length for Ongoing is a display
#   offset, not additional follow-up.

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
        patient = "patient",
        treatment = "treatment",
        start = "start",
        end = "end",
        response = "response",
        response_time = "response_time",
        stop = "stop"
    ),
    meta = c("line", "loh", "cdk", "endocrine", "platinum"),
    meta_titles = list(
        line = "Treatment line",
        loh = "BRCA2 zygosity",
        cdk = "CDK4/6i agent",
        endocrine = "Endocrine agent",
        platinum = "Previous platinum"
    ),
    palettes = list(
        treatment = "Qualitative.Dark2",
        response = "Quantitative.BluGrn",
        stop = "Qualitative.Safe",
        meta = "Qualitative.Prism"
    ),
    response_shown = c("CR", "PR", "SD"),
    treatment_order = c("CDK4/6i+ET", "Other", "PARPi"),
    sort = "PARPi",
    lane_half = 0.42,
    arrow_len = 2.2,
    meta_gap = 1.45,
    point_size = 5.5,
    stop_size = 6.0,
    meta_size = 7.2,
    labels = list(
        x = "Months",
        treatment = "Treatment",
        response = "Best response",
        stop = "Reason for stop"
    ),
    size = list(
        width = 11,
        height = 8.5
    )
)

load_packages(c("ggplot2", "ggprism", "ggnewscale", "patchwork", "readr", "grid"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
df <- as.data.frame(df, stringsAsFactors = FALSE)

p_col <- config$columns$patient
t_col <- config$columns$treatment
s_col <- config$columns$start
e_col <- config$columns$end
r_col <- config$columns$response
rt_col <- config$columns$response_time
st_col <- config$columns$stop

df[[p_col]] <- as.character(df[[p_col]])
df[[t_col]] <- as.character(df[[t_col]])
df[[s_col]] <- as.numeric(df[[s_col]])
df[[e_col]] <- as.numeric(df[[e_col]])
df[[r_col]] <- as.character(df[[r_col]])
df[[rt_col]] <- as.numeric(df[[rt_col]])
df[[st_col]] <- as.character(df[[st_col]])
if (any(!nzchar(df[[p_col]]) | !nzchar(df[[t_col]]))) {
    stop("patient and treatment must be non-empty.", call. = FALSE)
}
if (any(!is.finite(df[[s_col]]) | !is.finite(df[[e_col]]))) {
    stop("start and end must be finite.", call. = FALSE)
}
if (any(df[[e_col]] < df[[s_col]])) {
    stop("end must be >= start.", call. = FALSE)
}

meta_cols <- config$meta
meta_cols <- meta_cols[meta_cols %in% names(df)]
if (length(meta_cols)) {
    require_columns(df, meta_cols)
    for (nm in meta_cols) {
        df[[nm]] <- as.character(df[[nm]])
    }
}

patients <- unique(df[[p_col]])
if (nzchar(config$sort) && config$sort %in% df[[t_col]]) {
    sub <- df[df[[t_col]] == config$sort, , drop = FALSE]
    dur <- tapply(sub[[e_col]] - sub[[s_col]], sub[[p_col]], max)
    ord <- names(sort(dur, decreasing = TRUE))
    rest <- setdiff(patients, ord)
    patients <- c(ord, rest)
}
n_pt <- length(patients)
# First sorted patient at the top: last discrete level sits highest.
df[[p_col]] <- factor(df[[p_col]], levels = rev(patients))
df$y <- as.integer(df[[p_col]])

treat_lv <- unique(as.character(df[[t_col]]))
pref <- config$treatment_order
pref <- pref[pref %in% treat_lv]
treat_lv <- c(pref, setdiff(treat_lv, pref))
df[[t_col]] <- factor(df[[t_col]], levels = treat_lv)
n_treat <- length(treat_lv)
main_lv <- setdiff(treat_lv, "Other")
treat_cols <- palette_colors(config$palettes$treatment, n = max(length(main_lv), 1L))
names(treat_cols) <- main_lv
if ("Other" %in% treat_lv) {
    treat_cols <- c(
        treat_cols,
        Other = palette_colors("Qualitative.Safe")[[12]]
    )
}
treat_cols <- treat_cols[treat_lv]

resp_lv <- config$response_shown
ramp <- palette_colors(config$palettes$response)
resp_idx <- unique(round(seq(1, length(ramp), length.out = length(resp_lv))))
resp_cols <- ramp[resp_idx]
names(resp_cols) <- resp_lv

safe <- palette_colors(config$palettes$stop)
stop_cols <- c(
    Toxicity = safe[[length(safe)]],
    `Disease Progression` = safe[[5]],
    Ongoing = safe[[5]]
)

meta_map <- list()
if (length(meta_cols)) {
    pool <- palette_colors(config$palettes$meta)
    used <- character()
    for (nm in meta_cols) {
        lv <- unique(df[[nm]][nzchar(df[[nm]]) & !is.na(df[[nm]])])
        take <- pool[!.hex_key(pool) %in% .hex_key(used)]
        if (length(take) < length(lv)) {
            take <- palette_colors(config$palettes$meta, n = length(lv) + length(used))
            take <- take[!.hex_key(take) %in% .hex_key(used)]
        }
        cols <- take[seq_len(length(lv))]
        names(cols) <- lv
        meta_map[[nm]] <- cols
        used <- c(used, unname(cols))
    }
}

max_end <- tapply(df[[e_col]], df[[p_col]], max)
seg <- data.frame(
    y = as.integer(factor(names(max_end), levels = levels(df[[p_col]]))),
    xmax = as.numeric(max_end),
    stringsAsFactors = FALSE
)

resp <- df[
    !is.na(df[[r_col]]) &
        df[[r_col]] %in% resp_lv &
        is.finite(df[[rt_col]]),
    ,
    drop = FALSE
]
stop_pt <- df[
    !is.na(df[[st_col]]) &
        df[[st_col]] %in% c("Toxicity", "Disease Progression"),
    ,
    drop = FALSE
]
ongoing <- df[!is.na(df[[st_col]]) & df[[st_col]] == "Ongoing", , drop = FALSE]

n_meta <- length(meta_cols)
x_left <- if (n_meta) -config$meta_gap * (n_meta + 0.55) else 0
x_right <- max(df[[e_col]], na.rm = TRUE) + config$arrow_len + 0.8
half <- config$lane_half

meta_df <- NULL
hdr_df <- NULL
if (n_meta) {
    one <- df[!duplicated(df[[p_col]]), , drop = FALSE]
    blocks <- vector("list", n_meta)
    hdr <- vector("list", n_meta)
    for (i in seq_along(meta_cols)) {
        nm <- meta_cols[[i]]
        x0 <- -config$meta_gap * (n_meta - i + 1L)
        blocks[[i]] <- data.frame(
            x = x0,
            y = one$y,
            value = one[[nm]],
            category = nm,
            stringsAsFactors = FALSE
        )
        title <- config$meta_titles[[nm]]
        if (is.null(title) || !nzchar(title)) {
            title <- nm
        }
        hdr[[i]] <- data.frame(
            x = x0,
            y = 1 - half - 0.02,
            label = title,
            stringsAsFactors = FALSE
        )
    }
    meta_df <- do.call(rbind, blocks)
    hdr_df <- do.call(rbind, hdr)
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_segment(
        data = seg,
        aes(x = 0, xend = xmax, y = y, yend = y),
        colour = palette_colors("Qualitative.Safe")[[12]],
        linewidth = 0.35
    ) +
    geom_rect(
        data = df,
        aes(
            xmin = .data[[s_col]],
            xmax = .data[[e_col]],
            ymin = y - half,
            ymax = y + half,
            fill = .data[[t_col]]
        ),
        colour = NA
    ) +
    scale_fill_manual(
        name = config$labels$treatment,
        values = treat_cols,
        guide = guide_legend(order = 1)
    ) +
    new_scale_fill() +
    geom_point(
        data = resp,
        aes(x = .data[[rt_col]], y = y, fill = .data[[r_col]]),
        shape = 21,
        size = config$point_size,
        colour = palette_colors("Qualitative.Safe")[[5]],
        stroke = 0.35
    ) +
    scale_fill_manual(
        name = config$labels$response,
        values = resp_cols,
        breaks = resp_lv,
        guide = guide_legend(order = 2, override.aes = list(size = 5))
    ) +
    geom_point(
        data = stop_pt,
        aes(x = .data[[e_col]], y = y, colour = .data[[st_col]]),
        shape = 18,
        size = config$stop_size
    ) +
    geom_segment(
        data = ongoing,
        aes(
            x = .data[[e_col]],
            xend = .data[[e_col]] + config$arrow_len,
            y = y,
            yend = y,
            colour = .data[[st_col]]
        ),
        arrow = grid::arrow(length = unit(0.12, "cm"), type = "closed"),
        linewidth = 0.45
    ) +
    scale_colour_manual(
        name = config$labels$stop,
        values = stop_cols,
        guide = guide_legend(order = 3, override.aes = list(size = 5))
    )

if (!is.null(meta_df)) {
    for (i in seq_along(meta_cols)) {
        nm <- meta_cols[[i]]
        sub <- meta_df[meta_df$category == nm, , drop = FALSE]
        title <- config$meta_titles[[nm]]
        if (is.null(title) || !nzchar(title)) {
            title <- nm
        }
        p <- p +
            new_scale_fill() +
            geom_point(
                data = sub,
                aes(x = x, y = y, fill = value),
                shape = 22,
                size = config$meta_size,
                colour = palette_colors("Qualitative.Safe")[[5]],
                stroke = 0.35
            ) +
            scale_fill_manual(
                name = title,
                values = meta_map[[nm]],
                guide = "none"
            )
    }
    p <- p +
        geom_text(
            data = hdr_df,
            aes(x = x, y = y, label = label),
            size = 4.2,
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            colour = palette_colors("Qualitative.Safe")[[5]]
        )
}

ink <- palette_colors("Qualitative.Safe")[[5]]

p <- p +
    coord_cartesian(
        xlim = c(x_left, x_right),
        ylim = c(0.35, n_pt + 0.55),
        clip = "off"
    ) +
    scale_x_continuous(
        expand = c(0, 0),
        breaks = function(lims) {
            br <- pretty(c(0, max(lims[[2]], 0)), n = 5)
            br[br >= 0]
        },
        guide = ggplot2::guide_axis(cap = "both")
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = config$labels$x, y = NULL) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(3, 4, 32, 4, "mm"),
        axis.title.x = element_text(hjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none"
    )

theme_legend_box <- function(position) {
    theme_void() +
        theme(
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.position = position,
            legend.box = if (identical(position, "right")) "vertical" else "horizontal",
            legend.box.just = if (identical(position, "right")) "left" else "top",
            legend.title = element_text(size = 13, colour = ink),
            legend.text = element_text(size = 8.5, colour = ink),
            legend.key.size = unit(0.42, "cm"),
            legend.spacing.x = unit(0.45, "cm"),
            legend.spacing.y = unit(0.16, "cm"),
            legend.background = element_blank(),
            legend.box.background = element_blank(),
            legend.margin = margin(0, 2, 0, 2)
        )
}

p_event_leg <- ggplot() +
    geom_rect(
        data = df,
        aes(
            xmin = .data[[s_col]],
            xmax = .data[[e_col]],
            ymin = y - half,
            ymax = y + half,
            fill = .data[[t_col]]
        ),
        colour = NA
    ) +
    scale_fill_manual(
        name = config$labels$treatment,
        values = treat_cols,
        guide = guide_legend(order = 1, direction = "vertical")
    ) +
    new_scale_fill() +
    geom_point(
        data = resp,
        aes(x = .data[[rt_col]], y = y, fill = .data[[r_col]]),
        shape = 21,
        size = 3,
        colour = ink,
        stroke = 0.35
    ) +
    scale_fill_manual(
        name = config$labels$response,
        values = resp_cols,
        breaks = resp_lv,
        guide = guide_legend(
            order = 2,
            direction = "vertical",
            override.aes = list(size = 5)
        )
    ) +
    geom_point(
        data = stop_pt,
        aes(x = .data[[e_col]], y = y, colour = .data[[st_col]]),
        shape = 18,
        size = 3
    ) +
    geom_segment(
        data = ongoing,
        aes(
            x = .data[[e_col]],
            xend = .data[[e_col]] + config$arrow_len,
            y = y,
            yend = y,
            colour = .data[[st_col]]
        ),
        arrow = grid::arrow(length = unit(0.12, "cm"), type = "closed"),
        linewidth = 0.45
    ) +
    scale_colour_manual(
        name = config$labels$stop,
        values = stop_cols,
        guide = guide_legend(
            order = 3,
            direction = "vertical",
            override.aes = list(size = 5)
        )
    ) +
    theme_legend_box("right")

p_meta_leg <- ggplot()
if (!is.null(meta_df)) {
    for (i in seq_along(meta_cols)) {
        nm <- meta_cols[[i]]
        sub <- meta_df[meta_df$category == nm, , drop = FALSE]
        title <- config$meta_titles[[nm]]
        if (is.null(title) || !nzchar(title)) {
            title <- nm
        }
        if (i > 1L) {
            p_meta_leg <- p_meta_leg + new_scale_fill()
        }
        p_meta_leg <- p_meta_leg +
            geom_point(
                data = sub,
                aes(x = x, y = y, fill = value),
                shape = 22,
                size = 3,
                colour = ink,
                stroke = 0.35
            ) +
            scale_fill_manual(
                name = title,
                values = meta_map[[nm]],
                guide = guide_legend(
                    order = i,
                    direction = "vertical",
                    override.aes = list(size = 5.5, shape = 22)
                )
            )
    }
    p_meta_leg <- p_meta_leg + theme_legend_box("bottom")
}

guide_grob <- function(plt, prefer = "guide-box-bottom") {
    g <- ggplot2::ggplotGrob(plt)
    names <- g$layout$name
    prefer_hit <- which(names == prefer)
    hit <- unique(c(prefer_hit, grep("^guide-box", names)))
    for (i in hit) {
        grob <- g$grobs[[i]]
        if (!is.null(grob) && !inherits(grob, "zeroGrob")) {
            return(grob)
        }
    }
    grid::nullGrob()
}

tighten_guide <- function(grob) {
    if (!inherits(grob, "gtable")) {
        return(grob)
    }
    drop_null <- function(u) {
        n <- length(u)
        out <- grid::unit(rep(0, n), "cm")
        types <- grid::unitType(u)
        for (i in seq_len(n)) {
            if (!identical(types[[i]], "null")) {
                out[i] <- grid::convertUnit(u[i], "cm")
            }
        }
        out
    }
    grob$widths <- drop_null(grob$widths)
    grob$heights <- drop_null(grob$heights)
    grob
}

blank_bg <- theme(
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
)
wrap_leg <- function(grob) {
    patchwork::wrap_elements(full = grob) + blank_bg
}
anchor_grob <- function(grob, x, y, just) {
    grob <- tighten_guide(grob)
    grid::gTree(
        children = grid::gList(grob),
        vp = grid::viewport(
            x = grid::unit(x, "npc"),
            y = grid::unit(y, "npc"),
            just = just,
            width = grid::grobWidth(grob),
            height = grid::grobHeight(grob)
        )
    )
}

p <- p + patchwork::inset_element(
    wrap_leg(anchor_grob(
        guide_grob(p_event_leg, prefer = "guide-box-right"),
        x = 0.97,
        y = 0.04,
        just = c("right", "bottom")
    )),
    left = 0,
    bottom = 0,
    right = 1,
    top = 1,
    align_to = "panel",
    clip = FALSE,
    on_top = TRUE
)
if (!is.null(meta_df)) {
    tile_frac <- abs(x_left) / (x_right - x_left)
    p <- p + patchwork::inset_element(
        wrap_leg(anchor_grob(
            guide_grob(p_meta_leg, prefer = "guide-box-bottom"),
            x = 0.50,
            y = 1,
            just = c("centre", "top")
        )),
        left = tile_frac + 0.06,
        bottom = -0.28,
        right = 0.90,
        top = -0.12,
        align_to = "panel",
        clip = FALSE,
        on_top = TRUE
    )
}
p <- p + patchwork::plot_annotation(
    theme = theme(
        plot.background = element_rect(fill = NA, colour = NA),
        plot.margin = margin(0, 0, 0, 0, "mm")
    )
)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
suppressWarnings(save_ggplot(
    p, io$output,
    width = config$size$width,
    height = config$size$height
))
