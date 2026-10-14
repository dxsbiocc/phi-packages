#!/usr/bin/env Rscript
# Extract one PharmacoSet RDS into cells.csv, drugs.csv, profiles.csv, and
# raw.csv when a dose-response matrix is present. A plain list with those
# names is accepted so the Python tests can round-trip without PharmacoGx.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: Rscript extract_pset.R <input.rds> <output-dir>")
}
rds_path <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_df <- function(df, name) {
  if (is.null(df) || !nrow(df)) {
    return(invisible(NULL))
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  utils::write.csv(df, file.path(out_dir, name), row.names = FALSE, na = "")
}

pick_col <- function(df, candidates) {
  names_folded <- tolower(names(df))
  for (cand in candidates) {
    hit <- which(names_folded == tolower(cand))
    if (length(hit)) {
      return(df[[hit[[1]]]])
    }
  }
  NULL
}

norm_cells <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  ids <- rownames(df)
  source_model_id <- pick_col(df, c("source_model_id", "cellid", "sampleid"))
  if (is.null(source_model_id)) {
    source_model_id <- ids
  }
  cell_name <- pick_col(df, c("cell_name", "CellLine", "cellid", "sampleid"))
  if (is.null(cell_name)) {
    cell_name <- source_model_id
  }
  cellosaurus <- pick_col(df, c("cellosaurus", "Cellosaurus.Accession.id", "Cellosaurus.Accession"))
  tissue <- pick_col(df, c("tissue", "tissueid"))
  data.frame(
    source_model_id = as.character(source_model_id),
    cell_name = as.character(cell_name),
    cellosaurus = if (is.null(cellosaurus)) NA_character_ else as.character(cellosaurus),
    tissue = if (is.null(tissue)) NA_character_ else as.character(tissue),
    stringsAsFactors = FALSE
  )
}

norm_drugs <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  ids <- rownames(df)
  source_compound_id <- pick_col(df, c("source_compound_id", "drugid", "treatmentid"))
  if (is.null(source_compound_id)) {
    source_compound_id <- ids
  }
  name <- pick_col(df, c("name", "drug_name", "Drug.Name", "drugid"))
  cid <- pick_col(df, c("pubchem_cid", "cid", "PubChem.CID"))
  smiles <- pick_col(df, c("smiles", "SMILES", "smiles_cannonical"))
  inchikey <- pick_col(df, c("inchikey", "InChIKey"))
  brd <- pick_col(df, c("brd_id", "broad_id"))
  data.frame(
    source_compound_id = as.character(source_compound_id),
    name = if (is.null(name)) NA_character_ else as.character(name),
    pubchem_cid = if (is.null(cid)) NA_character_ else as.character(cid),
    pubchem_cid_alt = NA_character_,
    smiles = if (is.null(smiles)) NA_character_ else as.character(smiles),
    inchikey = if (is.null(inchikey)) NA_character_ else as.character(inchikey),
    brd_id = if (is.null(brd)) NA_character_ else as.character(brd),
    chembl_id = NA_character_,
    stringsAsFactors = FALSE
  )
}

melt_profiles <- function(info, profiles) {
  info <- as.data.frame(info, stringsAsFactors = FALSE)
  profiles <- as.data.frame(profiles, stringsAsFactors = FALSE)
  if (!("experiment_id" %in% names(info))) {
    info$experiment_id <- rownames(info)
  }
  if (!("experiment_id" %in% names(profiles))) {
    profiles$experiment_id <- rownames(profiles)
  }
  merged <- merge(info, profiles, by = "experiment_id", all = FALSE)
  cell_col <- NULL
  for (cand in c("source_model_id", "cellid", "sampleid")) {
    if (cand %in% names(merged)) {
      cell_col <- cand
      break
    }
  }
  drug_col <- NULL
  for (cand in c("source_compound_id", "drugid", "treatmentid")) {
    if (cand %in% names(merged)) {
      drug_col <- cand
      break
    }
  }
  if (is.null(cell_col) || is.null(drug_col)) {
    stop("sensitivityInfo has no cellid/drugid columns")
  }
  metric_names <- intersect(
    names(merged),
    c("auc", "ic50", "aac", "ec50", "auc_recomputed", "ic50_recomputed", "aac_recomputed", "ac50_recomputed", "ln_ic50")
  )
  if (!length(metric_names)) {
    stop("sensitivityProfiles has no auc/ic50/aac/ec50 columns")
  }
  pieces <- lapply(metric_names, function(metric) {
    data.frame(
      source_model_id = as.character(merged[[cell_col]]),
      source_compound_id = as.character(merged[[drug_col]]),
      metric = metric,
      value = suppressWarnings(as.numeric(merged[[metric]])),
      unit = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

obj <- readRDS(rds_path)

if (is.list(obj) && !methods::is(obj, "PharmacoSet") && !is.null(obj$cells)) {
  write_df(norm_cells(obj$cells), "cells.csv")
  if (!is.null(obj$drugs)) {
    write_df(norm_drugs(obj$drugs), "drugs.csv")
  }
  if (!is.null(obj$profiles)) {
    profiles <- as.data.frame(obj$profiles, stringsAsFactors = FALSE)
    if (all(c("source_model_id", "source_compound_id", "metric", "value") %in% names(profiles))) {
      write_df(profiles, "profiles.csv")
    } else if (!is.null(obj$sensitivity_info)) {
      write_df(melt_profiles(obj$sensitivity_info, profiles), "profiles.csv")
    }
  }
  if (!is.null(obj$raw)) {
    write_df(as.data.frame(obj$raw, stringsAsFactors = FALSE), "raw.csv")
  }
} else if (requireNamespace("PharmacoGx", quietly = TRUE) && methods::is(obj, "PharmacoSet")) {
  write_df(norm_cells(PharmacoGx::cellInfo(obj)), "cells.csv")
  write_df(norm_drugs(PharmacoGx::drugInfo(obj)), "drugs.csv")
  write_df(melt_profiles(PharmacoGx::sensitivityInfo(obj), PharmacoGx::sensitivityProfiles(obj)), "profiles.csv")
  raw <- tryCatch(PharmacoGx::sensitivityRaw(obj), error = function(err) NULL)
  if (is.data.frame(raw)) {
    write_df(raw, "raw.csv")
  } else {
    note <- "sensitivityRaw was not a data.frame; dose-response rows were not written."
    writeLines(note, file.path(out_dir, "extract_notes.txt"))
  }
} else {
  stop("RDS is not a list(cells, drugs, profiles) and PharmacoGx is not available to read a PharmacoSet. See the pipeline Dockerfile.")
}
