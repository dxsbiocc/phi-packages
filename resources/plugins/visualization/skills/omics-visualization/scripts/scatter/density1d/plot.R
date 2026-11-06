#!/usr/bin/env Rscript

# Template-ID: scatter-density1d
#
# Purpose:
#   Draw a 1D kernel density of a numeric value, optionally filled /
#   coloured by a group (identity, stack, or fill position), faceted,
#   or mirrored as a pair of series. style = "effect" draws a shared
#   reference density vs a panel-coloured case density with supplied
#   n / d / P annotations and optional mode lines.
#
# Inputs:
#   Default (example.tsv): carat / price / cut / depth.
#   Mirror (example_mirror.tsv): var1, var2.
#   Effect (example_effect.tsv): value, series, panel.
#   Optional annotations.tsv (effect): panel, n_label, stats.
#
# Output:
#   A PDF, PNG, or SVG 1D density figure.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for adjust, position, facet, and mirror columns.
#   Effect-size panels: style = "effect", series = reference|case,
#   panel colours from catalog Bold (blue / gold / pink). Do not
#   compute d or P here; annotation text is supplied. Histograms
#   stay scatter-hist. Ridgelines stay line-ridge. 2D densities
#   stay scatter-density2d.
#
# Scientific assumptions:
#   Each row is one observation. The kernel density is a display
#   smoothing of the supplied values, not a test. Mode lines mark
#   the plotted density peak only.

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
    # "single" | "mirror" | "effect"
    layout = "single",
    columns = list(
        x = "carat",
        group = NULL,            # e.g. "cut"
        series = "series",       # effect: "reference" | "case"
        facet = "panel",         # effect: panel / timepoint
        a = "var1",
        b = "var2"
    ),
    labels = list(
        title = NULL,
        x = "carat",
        y = NULL,
        a = "variable1",
        b = "variable2",
        caption = NULL
    ),
    adjust = 1,
    # "identity" | "stack" | "fill"
    position = "identity",
    # "density" | "count"  (count useful with position stack/fill)
    y_stat = "density",
    orientation = "x",           # "x" | "y"
    facet = FALSE,
    facet_scales = "fixed",      # "fixed" | "free_y" | "free"
    facet_nrow = NULL,
    facet_ncol = 1,
    alpha = 0.35,
    group_order = NULL,          # explicit level order for the group
    facet_order = NULL,
    linewidth = 0.5,
    xlim = NULL,                 # e.g. c(55, 70)
    # Effect layout: reference grey + per-panel case colours.
    effect = list(
        reference_level = "reference",
        case_level = "case",
        show_mode_lines = TRUE,
        mode_linewidth = 0.35,
        # Optional path beside the input, or absolute. Columns:
        # panel, n_label, stats (multi-line OK with \n).
        annotations = NULL,
        n_hjust = 0, n_vjust = 1.4, n_x = -Inf, n_y = Inf,
        stats_hjust = 1.05, stats_vjust = 1.4, stats_x = Inf, stats_y = Inf,
        label_size = 3.4
    ),
    palettes = list(
        group = "Qualitative.Bold",
        mirror_a = "Qualitative.Safe:7",
        mirror_b = "Qualitative.Safe:5",
        single = "Qualitative.Bold:3",
        # Soft grey for the shared diagnosis-free reference.
        reference = "Qualitative.Pastel / 3:4",
        # Case accents: Bold blue / gold / pink (indices 3,4,5).
        case = "Qualitative.Bold"
    ),
    theme = "bw",                # "bw" | "classic"
    base_size = 13,
    size = list(width = 7.5, height = 5)
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

order_levels <- function(vals, preferred) {
    lv <- unique(as.character(vals))
    if (is.null(preferred)) return(lv)
    c(intersect(as.character(preferred), lv), setdiff(lv, as.character(preferred)))
}

density_mode <- function(x, adjust = 1) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) return(NA_real_)
    d <- stats::density(x, adjust = adjust)
    d$x[[which.max(d$y)]]
}

df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
if (identical(config$layout, "mirror")) {
    require_columns(df, list(config$columns$a, config$columns$b))
    a <- as.numeric(df[[config$columns$a]])
    b <- as.numeric(df[[config$columns$b]])
    a <- a[is.finite(a)]
    b <- b[is.finite(b)]
    col_a <- pal_ref(config$palettes$mirror_a)
    col_b <- pal_ref(config$palettes$mirror_b)
    p <- ggplot() +
        geom_density(
            aes(x = a, y = after_stat(density)),
            fill = col_a, colour = col_a, alpha = 0.85,
            adjust = config$adjust, linewidth = config$linewidth
        ) +
        geom_density(
            aes(x = b, y = -after_stat(density)),
            fill = col_b, colour = col_b, alpha = 0.85,
            adjust = config$adjust, linewidth = config$linewidth
        ) +
        annotate(
            "label", x = Inf, y = Inf, label = config$labels$a,
            hjust = 1.1, vjust = 1.4, colour = col_a, size = 4,
            label.size = 0
        ) +
        annotate(
            "label", x = Inf, y = -Inf, label = config$labels$b,
            hjust = 1.1, vjust = -0.4, colour = col_b, size = 4,
            label.size = 0
        ) +
        geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
        labs(title = config$labels$title, x = config$labels$x, y = config$labels$y)
} else if (identical(config$layout, "effect")) {
    x_col <- config$columns$x
    s_col <- config$columns$series
    f_col <- config$columns$facet
    require_columns(df, list(x_col, s_col, f_col))
    df[[x_col]] <- as.numeric(df[[x_col]])
    df <- df[is.finite(df[[x_col]]), , drop = FALSE]
    if (!nrow(df)) stop("No finite values.", call. = FALSE)
    ec <- config$effect
    ref_lv <- as.character(ec$reference_level)
    case_lv <- as.character(ec$case_level)
    df[[s_col]] <- as.character(df[[s_col]])
    bad <- setdiff(unique(df[[s_col]]), c(ref_lv, case_lv))
    if (length(bad)) {
        stop(
            "effect series must be '", ref_lv, "' or '", case_lv,
            "'; got: ", paste(bad, collapse = ", "),
            call. = FALSE
        )
    }
    df[[f_col]] <- factor(df[[f_col]], levels = order_levels(df[[f_col]], config$facet_order))
    panels <- levels(df[[f_col]])
    case_cols <- expand_palette(palette_colors(config$palettes$case), length(panels))
    # Prefer Bold indices 3,4,5 (blue / gold / pink) when enough colours exist.
    bold <- palette_colors("Qualitative.Bold")
    if (identical(config$palettes$case, "Qualitative.Bold") && length(bold) >= 5L) {
        pick <- bold[c(3L, 4L, 5L, setdiff(seq_along(bold), c(3L, 4L, 5L)))]
        case_cols <- expand_palette(pick, length(panels))
    }
    names(case_cols) <- panels
    ref_col <- pal_ref(config$palettes$reference)

    # Colour key: reference + one colour per panel's case.
    fill_lv <- c(ref_lv, paste0(case_lv, "::", panels))
    fill_cols <- c(ref_col, unname(case_cols))
    names(fill_cols) <- fill_lv
    df$fill_key <- ifelse(
        df[[s_col]] == ref_lv,
        ref_lv,
        paste0(case_lv, "::", as.character(df[[f_col]]))
    )
    df$fill_key <- factor(df$fill_key, levels = fill_lv)

    p <- ggplot(df, aes(x = .data[[x_col]], fill = .data$fill_key, colour = .data$fill_key)) +
        geom_density(
            alpha = config$alpha,
            adjust = config$adjust,
            linewidth = config$linewidth,
            position = "identity"
        ) +
        scale_fill_manual(values = fill_cols, guide = "none") +
        scale_colour_manual(values = fill_cols, guide = "none") +
        facet_wrap(
            stats::as.formula(paste("~", f_col)),
            ncol = if (is.null(config$facet_ncol)) 1L else config$facet_ncol,
            nrow = config$facet_nrow,
            scales = config$facet_scales
        )

    if (isTRUE(ec$show_mode_lines)) {
        mode_rows <- list()
        for (pn in panels) {
            sub <- df[df[[f_col]] == pn, , drop = FALSE]
            for (ser in c(ref_lv, case_lv)) {
                xx <- sub[[x_col]][sub[[s_col]] == ser]
                mx <- density_mode(xx, adjust = config$adjust)
                if (!is.finite(mx)) next
                col <- if (identical(ser, ref_lv)) ref_col else case_cols[[pn]]
                mode_rows[[length(mode_rows) + 1L]] <- data.frame(
                    panel = factor(pn, levels = panels),
                    x = mx,
                    colour = col,
                    stringsAsFactors = FALSE
                )
            }
        }
        if (length(mode_rows)) {
            mode_df <- do.call(rbind, mode_rows)
            names(mode_df)[1] <- f_col
            for (i in seq_len(nrow(mode_df))) {
                p <- p + geom_vline(
                    data = mode_df[i, , drop = FALSE],
                    aes(xintercept = .data$x),
                    colour = mode_df$colour[[i]],
                    linewidth = ec$mode_linewidth,
                    inherit.aes = FALSE
                )
            }
        }
    }

    # Supplied annotations.
    ann_path <- ec$annotations
    if (!is.null(ann_path) && nzchar(as.character(ann_path))) {
        if (!grepl("^(/|[A-Za-z]:)", ann_path) && !file.exists(ann_path)) {
            cand <- file.path(dirname(io$input), ann_path)
            if (file.exists(cand)) ann_path <- cand
        }
        ann <- as.data.frame(read_table_auto(ann_path), stringsAsFactors = FALSE)
        require_columns(ann, list("panel"))
        for (nm in intersect(c("n_label", "stats"), names(ann))) {
            ann[[nm]] <- gsub("\\n", "\n", as.character(ann[[nm]]), fixed = TRUE)
        }
        ann$panel <- factor(as.character(ann$panel), levels = panels)
        if (any(is.na(ann$panel))) {
            stop("annotation panel levels must match facet levels: ",
                 paste(panels, collapse = ", "), call. = FALSE)
        }
        names(ann)[names(ann) == "panel"] <- f_col
        if ("n_label" %in% names(ann)) {
            ann$n_x <- ec$n_x
            ann$n_y <- ec$n_y
            for (i in seq_len(nrow(ann))) {
                pn <- as.character(ann[[f_col]][[i]])
                p <- p + geom_text(
                    data = ann[i, , drop = FALSE],
                    aes(x = .data$n_x, y = .data$n_y, label = .data$n_label),
                    hjust = ec$n_hjust, vjust = ec$n_vjust,
                    colour = case_cols[[pn]],
                    size = ec$label_size,
                    inherit.aes = FALSE
                )
            }
        }
        if ("stats" %in% names(ann)) {
            ann$stats_x <- ec$stats_x
            ann$stats_y <- ec$stats_y
            for (i in seq_len(nrow(ann))) {
                p <- p + geom_text(
                    data = ann[i, , drop = FALSE],
                    aes(x = .data$stats_x, y = .data$stats_y, label = .data$stats),
                    hjust = ec$stats_hjust, vjust = ec$stats_vjust,
                    colour = "grey20",
                    size = ec$label_size,
                    inherit.aes = FALSE
                )
            }
        }
    }

    if (!is.null(config$xlim)) {
        p <- p + coord_cartesian(xlim = config$xlim)
    }
    p <- p + labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        caption = config$labels$caption
    )
} else {
    x_col <- config$columns$x
    g_col <- config$columns$group
    require_columns(df, list(x_col))
    df[[x_col]] <- as.numeric(df[[x_col]])
    df <- df[is.finite(df[[x_col]]), , drop = FALSE]
    if (!nrow(df)) stop("No finite values.", call. = FALSE)
    use_group <- !is.null(g_col) && g_col %in% names(df)
    if (use_group) {
        df[[g_col]] <- factor(df[[g_col]], levels = order_levels(df[[g_col]], config$group_order))
        cols <- expand_palette(
            palette_colors(config$palettes$group), nlevels(df[[g_col]])
        )
        names(cols) <- levels(df[[g_col]])
    }
    count_y <- identical(config$y_stat, "count")
    horiz <- identical(config$orientation, "y")

    if (horiz) {
        p <- ggplot(df, aes(y = .data[[x_col]]))
        if (use_group) {
            mapping <- if (count_y) {
                aes(x = after_stat(count), fill = .data[[g_col]],
                    colour = .data[[g_col]])
            } else {
                aes(x = after_stat(density), fill = .data[[g_col]],
                    colour = .data[[g_col]])
            }
        } else if (count_y) {
            mapping <- aes(x = after_stat(count))
        } else {
            mapping <- aes(x = after_stat(density))
        }
    } else {
        p <- ggplot(df, aes(x = .data[[x_col]]))
        if (use_group) {
            mapping <- if (count_y) {
                aes(y = after_stat(count), fill = .data[[g_col]],
                    colour = .data[[g_col]])
            } else {
                aes(y = after_stat(density), fill = .data[[g_col]],
                    colour = .data[[g_col]])
            }
        } else if (count_y) {
            mapping <- aes(y = after_stat(count))
        } else {
            mapping <- aes(y = after_stat(density))
        }
    }

    dens_args <- list(
        mapping = mapping,
        adjust = config$adjust,
        alpha = config$alpha,
        linewidth = config$linewidth,
        position = if (use_group) config$position else "identity"
    )
    if (!use_group) {
        dens_args$fill <- pal_ref(config$palettes$single)
        dens_args$colour <- pal_ref(config$palettes$single)
    }
    p <- p + do.call(geom_density, dens_args)
    if (use_group) {
        p <- p +
            scale_fill_manual(values = cols, name = NULL) +
            scale_colour_manual(values = cols, name = NULL)
    }
    if (isTRUE(config$facet) && use_group) {
        p <- p + facet_wrap(
            stats::as.formula(paste("~", g_col)),
            scales = config$facet_scales
        )
    }
    if (!is.null(config$xlim) && !horiz) {
        p <- p + coord_cartesian(xlim = config$xlim)
    }
    p <- p + labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    )
}

if (identical(config$theme, "classic") || identical(config$layout, "effect")) {
    p <- p +
        theme_classic(base_size = config$base_size) +
        theme(
            axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            strip.background = element_blank(),
            strip.text = element_blank(),
            panel.spacing = grid::unit(0.85, "lines"),
            plot.caption = element_text(
                hjust = 0, size = rel(0.75), colour = "grey35",
                margin = margin(t = 8)
            ),
            legend.position = "none"
        )
} else {
    p <- p +
        theme_bw(base_size = config$base_size) +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
            axis.text = element_text(colour = "black"),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", colour = "black")
        )
}
plain_num <- function(b) format(b, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
if (identical(config$orientation, "y") && !identical(config$layout, "mirror") &&
    !identical(config$layout, "effect")) {
    p <- p + scale_x_continuous(labels = plain_num)
} else {
    p <- p + scale_y_continuous(labels = plain_num)
}
if (isTRUE(config$facet) && !identical(config$layout, "mirror") &&
    !identical(config$layout, "effect")) {
    p <- p + theme(legend.position = "none")
}
if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
    p <- p + theme(plot.title = element_blank())
} else if (identical(config$layout, "effect")) {
    p <- p + theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.05))
    )
}
if (is.null(config$labels$y) || !nzchar(as.character(config$labels$y))) {
    p <- p + theme(axis.title.y = element_blank())
}
if (is.null(config$labels$caption) || !nzchar(as.character(config$labels$caption))) {
    p <- p + theme(plot.caption = element_blank())
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width, height = config$size$height)
