# Cell-cell communication (CellPhoneDB-style heteromeric subunit-minimum permutation
# test on the unscaled RNA data layer). Expression co-occurrence indicates potential
# signaling capacity between cell types, NOT physical cell-cell contact.
need("Matrix", "ggplot2")

meta_cols <- colnames(obj[[]])
if (!(group_by %in% meta_cols)) {
  stop(sprintf("group column '%s' not found in metadata; available: %s",
               group_by, paste(meta_cols, collapse = ", ")), call. = FALSE)
}

groups_vec <- as.character(obj[[group_by]][, 1])
levels_vec <- sort(unique(groups_vec[!is.na(groups_vec) & nzchar(groups_vec)]))

if (length(levels_vec) < 2L) {
  info("communication skipped: fewer than 2 cell groups")
  analysis_notes <- c(analysis_notes, paste0("communication skipped: fewer than 2 groups in ", group_by))
  ccc_info <- list(skipped = "fewer than 2 groups", group_by = group_by)
} else {
  if (!is.null(interactions_path) && nzchar(interactions_path)) {
    if (!file.exists(interactions_path)) stop(paste("interactions_path not found:", interactions_path), call. = FALSE)
    sep_char <- if (grepl("[.]csv$", interactions_path, ignore.case = TRUE)) "," else "\t"
    lr_df <- utils::read.delim(interactions_path, sep = sep_char, stringsAsFactors = FALSE, check.names = FALSE)
    cols_low <- setNames(colnames(lr_df), tolower(trimws(colnames(lr_df))))
    if (all(c("ligand", "receptor") %in% names(cols_low))) {
      raw_pairs <- Map(list, as.character(lr_df[[cols_low[["ligand"]]]]), as.character(lr_df[[cols_low[["receptor"]]]]))
    } else if (all(c("source", "target") %in% names(cols_low))) {
      raw_pairs <- Map(list, as.character(lr_df[[cols_low[["source"]]]]), as.character(lr_df[[cols_low[["target"]]]]))
    } else {
      stop("interactions_path requires ligand/receptor or source/target columns", call. = FALSE)
    }
    resource_name <- basename(interactions_path)
  } else {
    raw_pairs <- list(
      c("TGFB1", "TGFBR1+TGFBR2"), c("TGFB2", "TGFBR1+TGFBR2"), c("BMP4", "BMPR1A+BMPR2"),
      c("VEGFA", "FLT1"), c("VEGFA", "KDR"), c("PDGFB", "PDGFRB"), c("EGF", "EGFR"),
      c("AREG", "EGFR"), c("HGF", "MET"), c("IGF1", "IGF1R"), c("FGF2", "FGFR1"),
      c("GAS6", "AXL"), c("MIF", "CD74+CXCR4"), c("MIF", "CD74+CD44"), c("APP", "CD74"),
      c("TNF", "TNFRSF1A"), c("TNF", "TNFRSF1B"), c("IL1B", "IL1R1+IL1RAP"),
      c("IL6", "IL6R+IL6ST"), c("IL10", "IL10RA+IL10RB"), c("IFNG", "IFNGR1+IFNGR2"),
      c("IL2", "IL2RA+IL2RB+IL2RG"), c("IL15", "IL15RA+IL2RB+IL2RG"), c("IL4", "IL4R"),
      c("IL7", "IL7R+IL2RG"), c("IL33", "IL1RL1"), c("CSF1", "CSF1R"), c("FLT3LG", "FLT3"),
      c("CCL2", "CCR2"), c("CCL3", "CCR1"), c("CCL4", "CCR5"), c("CCL5", "CCR5"),
      c("CXCL8", "CXCR1"), c("CXCL8", "CXCR2"), c("CXCL9", "CXCR3"), c("CXCL10", "CXCR3"),
      c("CXCL12", "CXCR4"), c("CXCL13", "CXCR5"), c("CCL19", "CCR7"), c("CCL21", "CCR7"),
      c("CX3CL1", "CX3CR1"), c("XCL1", "XCR1"),
      c("FN1", "ITGAV+ITGB1"), c("FN1", "ITGA5+ITGB1"), c("COL1A1", "ITGA1+ITGB1"),
      c("COL1A2", "ITGA2+ITGB1"), c("SPP1", "CD44"), c("SPP1", "ITGAV+ITGB1"),
      c("ICAM1", "ITGAL+ITGB2"), c("VCAM1", "ITGA4+ITGB1"), c("THBS1", "CD47"),
      c("LGALS9", "HAVCR2"), c("CD274", "PDCD1"), c("PDCD1LG2", "PDCD1"),
      c("CD80", "CD28"), c("CD86", "CD28"), c("CD80", "CTLA4"), c("CD86", "CTLA4"),
      c("CD40LG", "CD40"), c("FASLG", "FAS"), c("TNFSF10", "TNFRSF10B"),
      c("HLA-E", "KLRC1+KLRD1"), c("NECTIN2", "TIGIT"), c("PVR", "TIGIT"),
      c("DLL1", "NOTCH1"), c("JAG1", "NOTCH1"), c("JAG1", "NOTCH2"),
      c("WNT5A", "FZD5"), c("ANXA1", "FPR1"), c("C3", "C3AR1"), c("C5", "C5AR1")
    )
    resource_name <- "curated_starter_v1"
  }

  assay_name <- if (!is.null(assay) && nzchar(assay) && assay %in% SeuratObject::Assays(obj)) assay else SeuratObject::DefaultAssay(obj)
  expr_mat <- tryCatch(
    SeuratObject::LayerData(obj, assay = assay_name, layer = "data"),
    error = function(e) SeuratObject::GetAssayData(obj, assay = assay_name, slot = "data")
  )
  gene_names <- rownames(expr_mat)
  gene_lookup <- stats::setNames(gene_names, toupper(gene_names))

  resolved <- list()
  for (pr in raw_pairs) {
    l_subs <- trimws(unlist(strsplit(gsub("_", "+", pr[[1]]), "+", fixed = TRUE)))
    r_subs <- trimws(unlist(strsplit(gsub("_", "+", pr[[2]]), "+", fixed = TRUE)))
    l_subs <- l_subs[nzchar(l_subs)]
    r_subs <- r_subs[nzchar(r_subs)]
    if (length(l_subs) > 0L && length(r_subs) > 0L &&
        all(toupper(l_subs) %in% names(gene_lookup)) &&
        all(toupper(r_subs) %in% names(gene_lookup))) {
      l_matched <- unname(gene_lookup[toupper(l_subs)])
      r_matched <- unname(gene_lookup[toupper(r_subs)])
      resolved[[length(resolved) + 1L]] <- list(
        ligand = paste(l_matched, collapse = "+"),
        receptor = paste(r_matched, collapse = "+"),
        l_genes = l_matched,
        r_genes = r_matched
      )
    }
  }

  if (length(resolved) == 0L) {
    info("communication skipped: 0 ligand-receptor pairs matched features")
    analysis_notes <- c(analysis_notes, "communication skipped: no ligand-receptor subunits matched feature names")
    ccc_info <- list(skipped = "no ligand-receptor pairs matched", resource = resource_name, group_by = group_by)
  } else {
    uniq_genes <- sort(unique(unlist(lapply(resolved, function(x) c(x$l_genes, x$r_genes)))))
    sub_expr <- as.matrix(expr_mat[uniq_genes, , drop = FALSE])
    sub_pos <- (sub_expr > 0) * 1.0

    k_grp <- length(levels_vec)
    grp_idx <- match(groups_vec, levels_vec)
    valid_cells <- !is.na(grp_idx)
    sub_expr <- sub_expr[, valid_cells, drop = FALSE]
    sub_pos <- sub_pos[, valid_cells, drop = FALSE]
    grp_idx <- grp_idx[valid_cells]
    n_valid <- length(grp_idx)

    # Indicator matrix (cells x groups), column-normalized
    G_mat <- Matrix::sparseMatrix(i = seq_len(n_valid), j = grp_idx, x = 1.0, dims = c(n_valid, k_grp))
    grp_sizes <- as.numeric(Matrix::colSums(G_mat))
    grp_sizes[grp_sizes == 0] <- 1.0
    G_mean <- sweep(as.matrix(G_mat), 2L, grp_sizes, "/")

    obs_mean <- sub_expr %*% G_mean   # genes x k_grp
    obs_pct <- sub_pos %*% G_mean     # genes x k_grp

    # Filter pairs tested by subunit min expression pct
    tested_rows <- list()
    for (p_idx in seq_along(resolved)) {
      pr <- resolved[[p_idx]]
      l_m <- if (length(pr$l_genes) == 1L) obs_mean[pr$l_genes, ] else apply(obs_mean[pr$l_genes, , drop = FALSE], 2L, min)
      r_m <- if (length(pr$r_genes) == 1L) obs_mean[pr$r_genes, ] else apply(obs_mean[pr$r_genes, , drop = FALSE], 2L, min)
      l_p <- if (length(pr$l_genes) == 1L) obs_pct[pr$l_genes, ] else apply(obs_pct[pr$l_genes, , drop = FALSE], 2L, min)
      r_p <- if (length(pr$r_genes) == 1L) obs_pct[pr$r_genes, ] else apply(obs_pct[pr$r_genes, , drop = FALSE], 2L, min)
      for (si in seq_len(k_grp)) {
        if (l_p[[si]] < ccc_min_pct) next
        for (ri in seq_len(k_grp)) {
          if (r_p[[ri]] < ccc_min_pct) next
          score <- 0.5 * (l_m[[si]] + r_m[[ri]])
          if (score <= 0) next
          tested_rows[[length(tested_rows) + 1L]] <- list(
            p_idx = p_idx, si = si, ri = ri,
            sender = levels_vec[[si]], receiver = levels_vec[[ri]],
            ligand = pr$ligand, receptor = pr$receptor,
            ligand_mean = l_m[[si]], receptor_mean = r_m[[ri]],
            ligand_pct = l_p[[si]], receptor_pct = r_p[[ri]],
            score = score
          )
        }
      }
    }

    if (length(tested_rows) == 0L) {
      info("communication skipped: no sender-receiver pairs passed min_pct gate")
      analysis_notes <- c(analysis_notes, "communication skipped: no pairs passed min_pct")
      ccc_info <- list(skipped = "no pairs passed min_pct", resource = resource_name, group_by = group_by)
    } else {
      n_tested <- length(tested_rows)
      obs_scores <- vapply(tested_rows, function(x) x$score, numeric(1))
      exceed <- integer(n_tested)
      set.seed(seed)
      n_perms_i <- max(50L, as.integer(ccc_perms))
      for (perm in seq_len(n_perms_i)) {
        perm_idx <- sample.int(n_valid)
        pG <- G_mean[perm_idx, , drop = FALSE]
        p_mean <- sub_expr %*% pG
        p_scores <- numeric(n_tested)
        for (ti in seq_len(n_tested)) {
          tr <- tested_rows[[ti]]
          pr <- resolved[[tr$p_idx]]
          lm <- min(p_mean[pr$l_genes, tr$si])
          rm <- min(p_mean[pr$r_genes, tr$ri])
          p_scores[[ti]] <- 0.5 * (lm + rm)
        }
        exceed <- exceed + as.integer(p_scores >= obs_scores)
      }
      pvals <- (exceed + 1.0) / (n_perms_i + 1.0)
      fdrs <- stats::p.adjust(pvals, method = "BH")

      res_df <- data.frame(
        sender = vapply(tested_rows, `[[`, character(1), "sender"),
        receiver = vapply(tested_rows, `[[`, character(1), "receiver"),
        ligand = vapply(tested_rows, `[[`, character(1), "ligand"),
        receptor = vapply(tested_rows, `[[`, character(1), "receptor"),
        score = round(obs_scores, 5),
        ligand_mean = round(vapply(tested_rows, `[[`, numeric(1), "ligand_mean"), 5),
        receptor_mean = round(vapply(tested_rows, `[[`, numeric(1), "receptor_mean"), 5),
        ligand_pct = round(vapply(tested_rows, `[[`, numeric(1), "ligand_pct"), 4),
        receptor_pct = round(vapply(tested_rows, `[[`, numeric(1), "receptor_pct"), 4),
        pvalue = signif(pvals, 5),
        fdr = signif(fdrs, 5),
        significant = fdrs <= ccc_fdr,
        stringsAsFactors = FALSE
      )
      res_df <- res_df[order(res_df$fdr, -res_df$score), ]

      out_dir <- dirname(normalizePath(ccc_output, mustWork = FALSE))
      if (nzchar(out_dir) && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      utils::write.table(res_df, file = ccc_output, sep = "\t", quote = FALSE, row.names = FALSE)

      sig_df <- res_df[res_df$significant, , drop = FALSE]
      count_grid <- expand.grid(sender = levels_vec, receiver = levels_vec, stringsAsFactors = FALSE)
      count_grid$n_sig <- vapply(seq_len(nrow(count_grid)), function(i) {
        sum(sig_df$sender == count_grid$sender[[i]] & sig_df$receiver == count_grid$receiver[[i]])
      }, integer(1))
      p_heat <- ggplot2::ggplot(count_grid, ggplot2::aes(x = receiver, y = sender, fill = n_sig)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = n_sig), size = 3.2) +
        ggplot2::scale_fill_gradient(low = "#f8fafc", high = "#0284c7") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
        ggplot2::labs(title = "Significant L-R interactions", x = "Receiver", y = "Sender", fill = "Count")
      show_plot(p_heat, "communication-heatmap", width = 6.5, height = 5.5)

      top_plot <- utils::head(if (nrow(sig_df) > 0L) sig_df else res_df, 20L)
      top_plot$pair <- paste(top_plot$ligand, "->", top_plot$receptor)
      top_plot$edge <- paste(top_plot$sender, "->", top_plot$receiver)
      p_dot <- ggplot2::ggplot(top_plot, ggplot2::aes(x = edge, y = pair, size = score, color = -log10(pmax(fdr, 1e-6)))) +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient(low = "#94a3b8", high = "#e11d48") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 8)) +
        ggplot2::labs(title = "Top ligand-receptor interactions", x = "Sender -> Receiver", y = "Ligand -> Receptor", color = "-log10(FDR)", size = "Mean score")
      show_plot(p_dot, "communication-dotplot", width = 8.0, height = 5.5)

      analysis_steps <- c(analysis_steps, "communication")
      ccc_info <- list(
        resource = resource_name,
        group_by = group_by,
        n_groups = k_grp,
        n_pairs_tested = n_tested,
        n_significant = nrow(sig_df),
        fdr_cutoff = ccc_fdr,
        min_pct = ccc_min_pct,
        n_perms = n_perms_i,
        output_table = ccc_output,
        top_interactions = utils::head(res_df, 10L)
      )
      cat(sprintf("communication (%s): %d/%d significant (FDR <= %.3g)\n",
                  resource_name, nrow(sig_df), n_tested, ccc_fdr))
    }
  }
}
