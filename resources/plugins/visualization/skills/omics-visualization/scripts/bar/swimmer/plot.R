#!/usr/bin/env Rscript

# Template-ID: bar-swimmer
#
# Purpose:
#   Draw a clinical swimmer plot: one lane per patient, treatment
#   intervals as rectangles on a time axis, with response points,
#   stop marks, and optional left-hand annotation tiles.
#   style = "dose" draws dose-coloured lanes, black SD/PD glyphs,
#   escalate / compassionate marks, end labels, and an optional
#   x-axis break (white panel; preview_dose.png).
#
# Inputs:
#   One row per treatment interval. Default example:
#     - patient: lane id
#     - treatment: fill of the rectangle (CDK, Other, PARPi, ...)
#     - start, end: time in months (already aligned)
#     - response, response_time: best response and its scan time
#     - stop: Toxicity / Disease Progression / Ongoing
#     - line, loh, cdk, endocrine, platinum: optional annotation columns
#   Dose style (example_dose.tsv + events_dose.tsv):
#     - patient, dose, start, end, end_label
#     - events: patient, time, event (SD|PD|escalate|compassionate)
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
#   extra names in config$meta, not extra ids. Dose-lane symbols and
#   axis breaks are CONFIG (style = "dose"). Do not add an id per
#   Nature paper layout. Do not compute best response, progression,
#   or Kaplan-Meier curves. A two-point before/after display is
#   scatter-dumbbell. A population survival curve is line-survival.
#
# Scientific assumptions:
#   Start, end, response, and stop are supplied. This script does not
#   convert calendar dates, impute scans, or test treatments. Patient
#   order is a display sort. Arrow length for Ongoing is a display
#   offset, not additional follow-up. Axis break is a display device.

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
    # "lane": Nature-style multi-treatment lanes + left tiles (default).
    # "dose": dose-coloured bars, SD/PD glyphs, marks, end labels, axis break.
    style = "lane",
    columns = list(
        patient = "patient",
        treatment = "treatment",
        start = "start",
        end = "end",
        response = "response",
        response_time = "response_time",
        stop = "stop",
        dose = "dose",
        end_label = "end_label"
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
        meta = "Qualitative.Prism",
        # Dose fills: orange, purple, magenta, teal, gold, sky, grey(NE).
        dose = c(
            "Qualitative.Bold:7",
            "Qualitative.Bold:1",
            "Qualitative.Bold:5",
            "Qualitative.Bold:2",
            "Qualitative.Bold:4",
            "Qualitative.Safe:1",
            "Qualitative.Safe:12"
        )
    ),
    response_shown = c("CR", "PR", "SD"),
    treatment_order = c("CDK4/6i+ET", "Other", "PARPi"),
    dose_order = c(
        "DL1 (0.48 mg/kg)", "DL2 (1.44 mg/kg)", "DL3 (2.88 mg/kg)",
        "DL4 (4.32 mg/kg)", "DL5 (6.48 mg/kg)", "DL6 (9.72 mg/kg)", "NE"
    ),
    sort = "PARPi",
    # dose style: "duration" sorts longest lane at top.
    sort_dose = "duration",
    lane_half = 0.42,
    arrow_len = 2.2,
    meta_gap = 1.45,
    point_size = 5.5,
    stop_size = 6.0,
    meta_size = 7.2,
    dose = list(
        events = "events_dose.tsv",
        lane_half = 0.38,
        sd_size = 2.8,
        pd_size = 2.6,
        mark_size = 3.6,
        label_size = 3.0,
        label_nudge = 0.35,
        # Axis break: display gap between cut_low and cut_high.
        break_axis = TRUE,
        cut_low = 8,
        cut_high = 20,
        gap = 1.2,
        break_ticks = c(0, 2, 4, 6, 8, 20, 25)
    ),
    labels = list(
        x = "Months",
        treatment = "Treatment",
        response = "Best response",
        stop = "Reason for stop",
        dose_x = "Time from treatment initiation (months)",
        dose_y = "Tumor type and patient ID",
        dose_fill = "Dose Level",
        dose_symbols = "Response Symbols"
    ),
    size = list(
        width = 11,
        height = 8.5
    ),
    size_dose = list(
        width = 9.5,
        height = 8.2
    )
)

load_packages(c("ggplot2", "ggprism", "ggnewscale", "patchwork", "readr", "grid"))

pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols[[1]]
}

# -----------------------------------------------------------------------------
# DOSE STYLE
# -----------------------------------------------------------------------------
if (identical(config$style, "dose")) {
    df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
    p_col <- config$columns$patient
    d_col <- config$columns$dose
    s_col <- config$columns$start
    e_col <- config$columns$end
    lab_col <- config$columns$end_label
    require_columns(df, list(p_col, d_col, s_col, e_col))
    df[[p_col]] <- as.character(df[[p_col]])
    df[[d_col]] <- as.character(df[[d_col]])
    df[[s_col]] <- as.numeric(df[[s_col]])
    df[[e_col]] <- as.numeric(df[[e_col]])
    if (lab_col %in% names(df)) {
        df[[lab_col]] <- as.character(df[[lab_col]])
    } else {
        df[[lab_col]] <- ""
    }
    df <- df[is.finite(df[[s_col]]) & is.finite(df[[e_col]]), , drop = FALSE]
    if (!nrow(df)) stop("No finite dose lanes.", call. = FALSE)

    dc <- config$dose
    # Events
    ev_path <- dc$events
    if (!is.null(ev_path) && nzchar(as.character(ev_path))) {
        if (!grepl("^(/|[A-Za-z]:)", ev_path) && !file.exists(ev_path)) {
            cand <- file.path(dirname(io$input), ev_path)
            if (file.exists(cand)) ev_path <- cand
        }
        ev <- as.data.frame(read_table_auto(ev_path), stringsAsFactors = FALSE)
        require_columns(ev, list("patient", "time", "event"))
        ev$patient <- as.character(ev$patient)
        ev$time <- as.numeric(ev$time)
        ev$event <- as.character(ev$event)
        ev <- ev[is.finite(ev$time) & nzchar(ev$event), , drop = FALSE]
    } else {
        ev <- data.frame(
            patient = character(), time = numeric(), event = character(),
            stringsAsFactors = FALSE
        )
    }

    # Sort patients
    if (identical(config$sort_dose, "duration")) {
        dur <- tapply(df[[e_col]] - df[[s_col]], df[[p_col]], max)
        patients <- names(sort(dur, decreasing = TRUE))
        rest <- setdiff(unique(df[[p_col]]), patients)
        patients <- c(patients, rest)
    } else {
        patients <- unique(df[[p_col]])
    }
    df[[p_col]] <- factor(df[[p_col]], levels = rev(patients))
    df$y <- as.integer(df[[p_col]])

    dose_lv <- unique(as.character(df[[d_col]]))
    if (!is.null(config$dose_order)) {
        dose_lv <- c(
            intersect(as.character(config$dose_order), dose_lv),
            setdiff(dose_lv, as.character(config$dose_order))
        )
    } else {
        # Prefer DL# order, NE last.
        ne <- dose_lv[grepl("^NE$", dose_lv, ignore.case = TRUE)]
        other <- setdiff(dose_lv, ne)
        nums <- suppressWarnings(as.integer(sub("(?i)^DL\\s*(\\d+).*", "\\1", other, perl = TRUE)))
        other <- other[order(nums, other, na.last = TRUE)]
        dose_lv <- c(other, ne)
    }
    df[[d_col]] <- factor(df[[d_col]], levels = dose_lv)
    dose_cols <- vapply(seq_along(dose_lv), function(i) {
        refs <- config$palettes$dose
        pal_ref(refs[[((i - 1L) %% length(refs)) + 1L]])
    }, character(1))
    names(dose_cols) <- dose_lv

    # Axis compression for break
    compress <- function(x) {
        x <- as.numeric(x)
        if (!isTRUE(dc$break_axis)) return(x)
        out <- x
        hi <- x > dc$cut_low
        out[hi] <- dc$cut_low + dc$gap + (x[hi] - dc$cut_high)
        # Values in the gap (between cut_low and cut_high) squeeze to cut_low.
        mid <- x > dc$cut_low & x < dc$cut_high
        out[mid] <- dc$cut_low + dc$gap * (x[mid] - dc$cut_low) / (dc$cut_high - dc$cut_low)
        out
    }
    df$x0 <- compress(df[[s_col]])
    df$x1 <- compress(df[[e_col]])
    half <- dc$lane_half

    if (nrow(ev)) {
        ev <- ev[ev$patient %in% levels(df[[p_col]]), , drop = FALSE]
        ev$y <- as.integer(factor(ev$patient, levels = levels(df[[p_col]])))
        ev$x <- compress(ev$time)
        ev_sd <- ev[ev$event == "SD", , drop = FALSE]
        ev_pd <- ev[ev$event == "PD", , drop = FALSE]
        ev_esc <- ev[ev$event == "escalate", , drop = FALSE]
        ev_comp <- ev[ev$event == "compassionate", , drop = FALSE]
    } else {
        ev_sd <- ev_pd <- ev_esc <- ev_comp <- ev
    }

    lab_df <- df[nzchar(df[[lab_col]]), , drop = FALSE]
    lab_df$x <- lab_df$x1 + dc$label_nudge

    # Tick map
    ticks <- dc$break_ticks
    tick_pos <- compress(ticks)
    x_max <- max(c(df$x1, tick_pos), na.rm = TRUE) + dc$label_nudge + 1.2

    p <- ggplot() +
        geom_rect(
            data = df,
            aes(
                xmin = .data$x0, xmax = .data$x1,
                ymin = .data$y - half, ymax = .data$y + half,
                fill = .data[[d_col]]
            ),
            colour = NA
        ) +
        scale_fill_manual(
            name = config$labels$dose_fill,
            values = dose_cols,
            guide = guide_legend(order = 2, override.aes = list(colour = NA))
        )

    if (nrow(ev_sd)) {
        p <- p + geom_point(
            data = ev_sd,
            aes(x = .data$x, y = .data$y),
            shape = 16, size = dc$sd_size, colour = "grey10"
        )
    }
    if (nrow(ev_pd)) {
        p <- p + geom_point(
            data = ev_pd,
            aes(x = .data$x, y = .data$y),
            shape = 15, size = dc$pd_size, colour = "grey10"
        )
    }
    if (nrow(ev_esc)) {
        p <- p + geom_text(
            data = ev_esc,
            aes(x = .data$x, y = .data$y + half * 0.55),
            label = "*", size = dc$mark_size, colour = "grey10",
            vjust = 0, fontface = "bold"
        )
    }
    if (nrow(ev_comp)) {
        p <- p + geom_text(
            data = ev_comp,
            aes(x = .data$x, y = .data$y + half * 0.55),
            label = "#", size = dc$mark_size, colour = "grey10",
            vjust = 0, fontface = "bold"
        )
    }
    if (nrow(lab_df)) {
        p <- p + geom_text(
            data = lab_df,
            aes(x = .data$x, y = .data$y, label = .data[[lab_col]]),
            hjust = 0, size = dc$label_size, colour = "grey15"
        )
    }

    # Symbol legend via invisible dummy points (shape codes for * / #).
    sym_df <- data.frame(
        x = NA_real_, y = NA_real_,
        symbol = factor(
            c(
                "Stable Disease (SD)",
                "Progressive Disease (PD)",
                "Escalated to next dose level",
                "Compassionate use"
            ),
            levels = c(
                "Stable Disease (SD)",
                "Progressive Disease (PD)",
                "Escalated to next dose level",
                "Compassionate use"
            )
        ),
        stringsAsFactors = FALSE
    )
    p <- p +
        geom_point(
            data = sym_df,
            aes(x = .data$x, y = .data$y, shape = .data$symbol),
            na.rm = TRUE, size = 3, colour = "grey10"
        ) +
        scale_shape_manual(
            name = config$labels$dose_symbols,
            values = c(
                "Stable Disease (SD)" = 16,
                "Progressive Disease (PD)" = 15,
                "Escalated to next dose level" = 42,  # *
                "Compassionate use" = 35              # #
            ),
            guide = guide_legend(order = 1, override.aes = list(size = 3.8))
        )

    p <- p +
        scale_y_continuous(
            breaks = df$y[!duplicated(df[[p_col]])],
            labels = as.character(df[[p_col]][!duplicated(df[[p_col]])]),
            expand = expansion(add = 0.6)
        ) +
        scale_x_continuous(
            breaks = tick_pos,
            labels = as.character(ticks),
            limits = c(0, x_max),
            expand = expansion(mult = c(0, 0.02))
        ) +
        labs(
            x = config$labels$dose_x,
            y = config$labels$dose_y
        ) +
        theme_classic(base_size = 12) +
        theme(
            axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.position = "inside",
            legend.position.inside = c(0.98, 0.02),
            legend.justification = c(1, 0),
            legend.box = "vertical",
            legend.background = element_rect(fill = "white", colour = NA),
            legend.key = element_rect(fill = NA, colour = NA),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.margin = margin(8, 10, 8, 8)
        )

    if (isTRUE(dc$break_axis)) {
        # Small break marks on the x axis near the gap.
        br_x <- dc$cut_low + dc$gap / 2
        p <- p +
            annotate(
                "text", x = br_x, y = -0.15,
                label = "//", size = 4, colour = "grey30"
            ) +
            coord_cartesian(ylim = c(0.4, max(df$y) + 0.6), clip = "off")
    }

    save_ggplot(
        p, io$output,
        width = config$size_dose$width,
        height = config$size_dose$height
    )
    quit(save = "no")
}

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
df <- as.data.frame(df, stringsAsFactors = FALSE)

p_col <- config$columns$patient
t_col <- config$columns$treatment
s_col <- config$columns$start
e_col <- config$columns$end
r_col <- config$columns$response
rt_col <- config$columns$response_time
st_col <- config$columns$stop
require_columns(df, list(p_col, t_col, s_col, e_col, r_col, rt_col, st_col))

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
