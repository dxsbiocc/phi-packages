#!/usr/bin/env Rscript

# Layout-ID: wide-top-two-bottom
# Purpose: place a wide overview above two complementary detail panels.
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
p_overview <- demo$heatmap
p_detail_a <- demo$scatter
p_detail_b <- demo$bar

# Demo previews hide legends so the layout geometry remains visible. For real
# panels, remove this block and choose an intentional keep/collect strategy.
panel_theme <- ggplot2::theme(legend.position = "none")
p_overview <- p_overview + panel_theme
p_detail_a <- p_detail_a + panel_theme
p_detail_b <- p_detail_b + panel_theme

# --- COMPOSITION: edit the visible layout code directly ----------------------

design <- "
AAAA
BBCC
"

figure <- patchwork::wrap_plots(
    A = p_overview,
    B = p_detail_a,
    C = p_detail_b,
    design = design
) +
    patchwork::plot_layout(guides = "keep", heights = c(1.15, 1)) +
    patchwork::plot_annotation(tag_levels = "a")

width <- 8.2
height <- 6
audit_path <- layout_audit_path(output_path)

audit_patchwork_layout(
    figure,
    layout_id = "wide-top-two-bottom",
    json_path = audit_path,
    width = width,
    height = height
)
save_layout_figure(figure, output_path, width = width, height = height)
message("Audit: ", normalizePath(audit_path, mustWork = TRUE))
