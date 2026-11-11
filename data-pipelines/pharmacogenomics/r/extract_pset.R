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
  # Experiments are keyed by row name. An experiment_id column in sensitivityInfo is not that
  # key (CTRPv2 uses it for the replicate number), so it is never used to join. When the row
  # names do not line up but both tables have the same rows, they are joined by position.
  if (!length(intersect(rownames(info), rownames(profiles))) && nrow(info) == nrow(profiles)) {
    rownames(profiles) <- rownames(info)
  }
  info[[".pgx_key"]] <- rownames(info)
  profiles[[".pgx_key"]] <- rownames(profiles)
  profiles <- profiles[, setdiff(names(profiles), names(info)[names(info) != ".pgx_key"]), drop = FALSE]
  merged <- merge(info, profiles, by = ".pgx_key", all = FALSE)
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
  wanted <- c("auc", "ic50", "aac", "ec50", "auc_recomputed", "ic50_recomputed", "aac_recomputed", "ac50_recomputed", "ln_ic50")
  metric_names <- names(merged)[tolower(names(merged)) %in% wanted]
  if (!length(metric_names)) {
    stop("sensitivityProfiles has no auc/ic50/aac/ec50 columns")
  }
  pieces <- lapply(metric_names, function(metric) {
    data.frame(
      source_model_id = as.character(merged[[cell_col]]),
      source_compound_id = as.character(merged[[drug_col]]),
      metric = tolower(metric),
      value = suppressWarnings(as.numeric(merged[[metric]])),
      unit = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

# A PharmacoSet saved before PharmacoGx 3 keeps its tables in the slots cell, drug and
# sensitivity, and PharmacoGx 3 no longer exports cellInfo() and drugInfo() for it. Reading
# those slots directly works with either version and does not need PharmacoGx at all. A
# PharmacoSet saved by PharmacoGx 3 uses sampleInfo(), treatmentInfo() and friends.
pset_tables <- function(obj) {
  sensitivity <- attr(obj, "sensitivity")
  if (!is.null(sensitivity) && !is.null(attr(obj, "cell")) && !is.null(attr(obj, "drug"))) {
    return(list(
      cells = attr(obj, "cell"),
      drugs = attr(obj, "drug"),
      info = sensitivity$info,
      profiles = sensitivity$profiles,
      raw = sensitivity$raw
    ))
  }
  if (!requireNamespace("PharmacoGx", quietly = TRUE)) {
    stop("This PharmacoSet needs PharmacoGx to be read. See the pipeline Dockerfile.")
  }
  list(
    cells = PharmacoGx::sampleInfo(obj),
    drugs = PharmacoGx::treatmentInfo(obj),
    info = PharmacoGx::sensitivityInfo(obj),
    profiles = PharmacoGx::sensitivityProfiles(obj),
    raw = tryCatch(PharmacoGx::sensitivityRaw(obj), error = function(err) NULL)
  )
}

# sensitivityRaw is an array of experiment x dose slot x (Dose, Viability). Melt it to one
# row per measured dose. Set PGX_PSET_RAW=0 to skip it for very large PharmacoSets.
melt_raw <- function(raw, info) {
  if (identical(Sys.getenv("PGX_PSET_RAW"), "0")) {
    return(NULL)
  }
  if (is.data.frame(raw)) {
    return(raw)
  }
  if (!is.array(raw) || length(dim(raw)) != 3 || dim(raw)[[3]] < 2) {
    return(NULL)
  }
  info <- as.data.frame(info, stringsAsFactors = FALSE)
  ids <- dimnames(raw)[[1]]
  if (is.null(ids)) {
    ids <- rownames(info)
  }
  cell <- NULL
  drug <- NULL
  for (cand in c("source_model_id", "cellid", "sampleid")) {
    if (cand %in% names(info)) { cell <- info[[cand]]; break }
  }
  for (cand in c("source_compound_id", "drugid", "treatmentid")) {
    if (cand %in% names(info)) { drug <- info[[cand]]; break }
  }
  at <- match(ids, rownames(info))
  if (is.null(cell) || is.null(drug) || anyNA(at)) {
    return(NULL)
  }
  dose <- suppressWarnings(matrix(as.numeric(raw[, , 1]), nrow = dim(raw)[[1]]))
  viability <- suppressWarnings(matrix(as.numeric(raw[, , 2]), nrow = dim(raw)[[1]]))
  hit <- which(!is.na(dose) & !is.na(viability), arr.ind = TRUE)
  if (!nrow(hit)) {
    return(NULL)
  }
  data.frame(
    source_model_id = as.character(cell[at[hit[, 1]]]),
    source_compound_id = as.character(drug[at[hit[, 1]]]),
    dose = dose[hit],
    value = viability[hit],
    metric = "VIABILITY",
    unit = "percent",
    stringsAsFactors = FALSE
  )
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
} else if (inherits(obj, "PharmacoSet")) {
  tables <- pset_tables(obj)
  write_df(norm_cells(tables$cells), "cells.csv")
  write_df(norm_drugs(tables$drugs), "drugs.csv")
  write_df(melt_profiles(tables$info, tables$profiles), "profiles.csv")
  raw <- melt_raw(tables$raw, tables$info)
  if (is.data.frame(raw)) {
    write_df(raw, "raw.csv")
  } else {
    note <- "sensitivityRaw was not usable; dose-response rows were not written."
    writeLines(note, file.path(out_dir, "extract_notes.txt"))
  }
} else {
  stop("RDS is neither a PharmacoSet nor a list(cells, drugs, profiles).")
}
