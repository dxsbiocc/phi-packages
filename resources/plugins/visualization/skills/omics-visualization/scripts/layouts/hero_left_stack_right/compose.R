#!/usr/bin/env Rscript

# Layout-ID: hero-left-stack-right
# Purpose: place decisive evidence on the left and two supporting panels right.
# Usage: Rscript compose.R <output.{pdf,png,svg}>
# Adaptation: replace PANEL OBJECTS, then edit COMPOSITION directly.

local({
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (!length(file_arg)) stop("Run this file with Rscript.", call. = FALSE)
    dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
    for (i in seq_len(8)) {
        candidate <- file.path(dir, "lib", "layout_common.R")
        if (file.exists(candidate)) {
            source(normalizePath(candidate), chdir = FALSE)
            return(invisible())
        }
        parent <- dirname(dir)
        if (identical(parent, dir)) break
        dir <- parent
    }
    stop("Cannot find scripts/layouts/lib/layout_common.R", call. = FALSE)
})

load_layout_packages()
output_path <- parse_layout_output()

# --- PANEL OBJECTS: replace with project-local plot objects ------------------

demo <- layout_demo_panels()
p_main <- demo$scatter
p_control <- demo$bar
p_validation <- demo$line

# Demo previews hide legends so the layout geometry remains visible. For real
# panels, remove this block and choose an intentional keep/collect strategy.
panel_theme <- ggplot2::theme(legend.position = "none")
p_main <- p_main + panel_theme
p_control <- p_control + panel_theme
p_validation <- p_validation + panel_theme

# --- COMPOSITION: edit the visible layout code directly ----------------------

design <- "
AAB
AAC
"

figure <- patchwork::wrap_plots(
    A = p_main,
    B = p_control,
    C = p_validation,
    design = design
) +
    patchwork::plot_layout(guides = "keep") +
    patchwork::plot_annotation(tag_levels = "a")

width <- 8.2
height <- 5.4
audit_path <- layout_audit_path(output_path)

audit_patchwork_layout(
    figure,
    layout_id = "hero-left-stack-right",
    json_path = audit_path,
    width = width,
    height = height
)
save_layout_figure(figure, output_path, width = width, height = height)
message("Audit: ", normalizePath(audit_path, mustWork = TRUE))
