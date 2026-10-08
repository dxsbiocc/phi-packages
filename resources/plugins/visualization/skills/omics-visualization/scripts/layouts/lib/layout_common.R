# Shared execution and rendered-geometry helpers for canonical layout templates.

load_layout_packages <- function() {
    packages <- c("ggplot2", "patchwork", "jsonlite")
    missing <- packages[!vapply(
        packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )]
    if (length(missing) > 0L) {
        stop(
            paste(
                "Missing R package(s):", paste(missing, collapse = ", "),
                "\nInstall them before running this layout."
            ),
            call. = FALSE
        )
    }
    suppressPackageStartupMessages({
        for (package in packages) {
            library(package, character.only = TRUE)
        }
    })
    invisible(packages)
}

parse_layout_output <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    script <- if (length(script_arg)) {
        basename(sub("^--file=", "", script_arg[[1]]))
    } else {
        "compose.R"
    }
    usage <- sprintf("Usage: Rscript %s <output.{pdf,png,svg}>", script)

    if (length(args) == 1L && args[[1]] %in% c("-h", "--help")) {
        cat(usage, "\n", sep = "")
        quit(save = "no", status = 0L)
    }
    if (length(args) != 1L) {
        stop(usage, call. = FALSE)
    }
    args[[1]]
}

layout_audit_path <- function(output_path) {
    extension <- tools::file_ext(output_path)
    if (!nzchar(extension)) {
        return(paste0(output_path, ".layout-audit.json"))
    }
    sub(
        paste0("\\.", extension, "$"),
        ".layout-audit.json",
        output_path,
        ignore.case = TRUE
    )
}

layout_demo_panels <- function() {
    set.seed(2026L)
    palette <- c(control = "#5F8072", treated = "#C9785E")
    base_theme <- ggplot2::theme_minimal(base_size = 9) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold", size = 9),
            legend.position = "bottom",
            legend.title = ggplot2::element_blank(),
            plot.background = ggplot2::element_rect(
                fill = "white",
                colour = NA
            )
        )

    scatter_data <- data.frame(
        x = c(rnorm(18, -0.4, 0.65), rnorm(18, 0.65, 0.6)),
        y = c(rnorm(18, 0.2, 0.55), rnorm(18, 0.9, 0.55)),
        group = rep(c("control", "treated"), each = 18)
    )
    p_scatter <- ggplot2::ggplot(
        scatter_data,
        ggplot2::aes(x = x, y = y, colour = group)
    ) +
        ggplot2::geom_point(size = 1.8, alpha = 0.82) +
        ggplot2::scale_colour_manual(values = palette) +
        ggplot2::labs(title = "Primary pattern", x = "Component 1", y = "Component 2") +
        base_theme

    bar_data <- data.frame(
        pathway = factor(c("Repair", "Cycle", "Signal", "Metabolism"),
            levels = rev(c("Repair", "Cycle", "Signal", "Metabolism"))
        ),
        score = c(3.8, 3.1, 2.5, 1.8)
    )
    p_bar <- ggplot2::ggplot(
        bar_data,
        ggplot2::aes(x = score, y = pathway)
    ) +
        ggplot2::geom_col(width = 0.68, fill = "#D4A84F") +
        ggplot2::labs(title = "Decomposition", x = "Score", y = NULL) +
        base_theme +
        ggplot2::theme(legend.position = "none")

    line_data <- expand.grid(
        time = 1:5,
        group = c("control", "treated")
    )
    line_data$value <- c(1.0, 1.3, 1.5, 1.65, 1.8, 1.0, 1.5, 2.0, 2.5, 2.8)
    p_line <- ggplot2::ggplot(
        line_data,
        ggplot2::aes(x = time, y = value, colour = group)
    ) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 1.5) +
        ggplot2::scale_colour_manual(values = palette) +
        ggplot2::labs(title = "Validation", x = "Time", y = "Response") +
        base_theme

    heatmap_data <- expand.grid(row = LETTERS[1:5], column = paste0("S", 1:5))
    heatmap_data$value <- sin(seq_len(nrow(heatmap_data)) / 3)
    p_heatmap <- ggplot2::ggplot(
        heatmap_data,
        ggplot2::aes(x = column, y = row, fill = value)
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
            low = "#668DA1",
            mid = "#F5F0E4",
            high = "#B8674B",
            midpoint = 0
        ) +
        ggplot2::labs(title = "Boundary", x = NULL, y = NULL, fill = "Value") +
        base_theme +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.position = "right"
        )

    list(
        scatter = p_scatter,
        bar = p_bar,
        line = p_line,
        heatmap = p_heatmap
    )
}

measure_patchwork_panels <- function(figure, width, height) {
    grDevices::pdf(NULL, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)

    grid::grid.newpage()
    grid::grid.draw(patchwork::patchworkGrob(figure))
    grid::grid.force()

    viewport_names <- grid::grid.ls(
        viewports = TRUE,
        grobs = FALSE,
        print = FALSE
    )$name
    viewport_names <- unique(viewport_names[
        grepl("^panel-[0-9]+\\.", viewport_names)
    ])
    if (!length(viewport_names)) {
        return(data.frame())
    }

    measured <- lapply(viewport_names, function(viewport_name) {
        grid::seekViewport(viewport_name)
        lower <- grid::deviceLoc(
            grid::unit(0, "npc"),
            grid::unit(0, "npc"),
            valueOnly = TRUE
        )
        upper <- grid::deviceLoc(
            grid::unit(1, "npc"),
            grid::unit(1, "npc"),
            valueOnly = TRUE
        )
        grid::upViewport(0)

        data.frame(
            panel = as.integer(sub("^panel-([0-9]+).*", "\\1", viewport_name)),
            left = lower$x * 72,
            right = upper$x * 72,
            bottom = lower$y * 72,
            top = upper$y * 72,
            width = (upper$x - lower$x) * 72,
            height = (upper$y - lower$y) * 72,
            stringsAsFactors = FALSE
        )
    })

    result <- do.call(rbind, measured)
    result[order(result$panel), , drop = FALSE]
}

audit_patchwork_layout <- function(
    figure,
    layout_id,
    json_path,
    width,
    height,
    tolerance_pt = 1.5,
    strict = TRUE
) {
    panels <- measure_patchwork_panels(figure, width, height)
    expected <- switch(
        layout_id,
        "equal-2x2" = 4L,
        "hero-left-stack-right" = 3L,
        "wide-top-two-bottom" = 3L,
        stop("Unknown layout id: ", layout_id, call. = FALSE)
    )

    checks <- list()
    add_spread <- function(name, values) {
        delta <- if (length(values)) diff(range(values)) else Inf
        checks[[length(checks) + 1L]] <<- list(
            name = name,
            delta_pt = unname(delta),
            passed = is.finite(delta) && delta <= tolerance_pt
        )
    }
    add_delta <- function(name, first, second) {
        delta <- abs(first - second)
        checks[[length(checks) + 1L]] <<- list(
            name = name,
            delta_pt = unname(delta),
            passed = is.finite(delta) && delta <= tolerance_pt
        )
    }

    if (nrow(panels) == expected) {
        panel <- function(index) panels[panels$panel == index, , drop = FALSE]

        if (layout_id == "equal-2x2") {
            add_spread("peer panel widths", panels$width)
            add_spread("peer panel heights", panels$height)
            add_delta("left column alignment", panel(1)$left, panel(3)$left)
            add_delta("right column alignment", panel(2)$right, panel(4)$right)
            add_delta("top row alignment", panel(1)$top, panel(2)$top)
            add_delta("bottom row alignment", panel(3)$bottom, panel(4)$bottom)
            add_delta(
                "horizontal gutter consistency",
                panel(2)$left - panel(1)$right,
                panel(4)$left - panel(3)$right
            )
            add_delta(
                "vertical gutter consistency",
                panel(1)$bottom - panel(3)$top,
                panel(2)$bottom - panel(4)$top
            )
        } else if (layout_id == "hero-left-stack-right") {
            add_delta("supporting panel widths", panel(2)$width, panel(3)$width)
            add_delta("supporting panel heights", panel(2)$height, panel(3)$height)
            add_delta("supporting left edges", panel(2)$left, panel(3)$left)
            add_delta("supporting right edges", panel(2)$right, panel(3)$right)
            add_delta("hero and upper panel top", panel(1)$top, panel(2)$top)
            add_delta("hero and lower panel bottom", panel(1)$bottom, panel(3)$bottom)
        } else if (layout_id == "wide-top-two-bottom") {
            add_delta("detail panel widths", panel(2)$width, panel(3)$width)
            add_delta("detail panel heights", panel(2)$height, panel(3)$height)
            add_delta("detail panel tops", panel(2)$top, panel(3)$top)
            add_delta("detail panel bottoms", panel(2)$bottom, panel(3)$bottom)
            add_delta("overview and lower-left edge", panel(1)$left, panel(2)$left)
            add_delta("overview and lower-right edge", panel(1)$right, panel(3)$right)
        }
    }

    auditable <- nrow(panels) == expected && length(checks) > 0L
    passed <- auditable && all(vapply(checks, `[[`, logical(1), "passed"))
    status <- if (!auditable) {
        "NOT_AUDITABLE"
    } else if (passed) {
        "PASS"
    } else {
        "FAIL"
    }

    report <- list(
        schema_version = "omics-visualization.layout-audit.v1",
        layout_id = layout_id,
        status = status,
        tolerance_pt = tolerance_pt,
        output_width_in = width,
        output_height_in = height,
        panels = panels,
        checks = checks
    )

    dir.create(dirname(json_path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(
        report,
        json_path,
        auto_unbox = TRUE,
        pretty = TRUE,
        dataframe = "rows",
        digits = 5
    )

    if (strict && status != "PASS") {
        stop(
            sprintf("Layout audit %s; see %s", status, json_path),
            call. = FALSE
        )
    }
    invisible(report)
}

save_layout_figure <- function(figure, path, width, height, dpi = 300) {
    format <- tolower(tools::file_ext(path))
    if (!format %in% c("pdf", "png", "svg")) {
        stop("Unsupported output format: ", format, call. = FALSE)
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
        filename = path,
        plot = figure,
        width = width,
        height = height,
        units = "in",
        dpi = dpi,
        bg = "transparent"
    )
    message("Saved: ", normalizePath(path, mustWork = TRUE))
    invisible(path)
}
