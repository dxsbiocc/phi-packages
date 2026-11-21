#!/usr/bin/env Rscript
# Report which ChIP/CUT&Tag downstream packages this environment provides, and
# which analysis presets can therefore run here. Phi environments are read-only
# and prebuilt: a missing package is a fact to report with a concrete next step
# (project environment via env_request), never an in-environment install.
#
# Presets mirror assets/steps.json: 1 TF/narrow DiffBind+DESeq2, 2 broad csaw,
# 3 CUT&Tag DiffBind+edgeR, 4 consensus, 5 annotation, 6 count exploration.

args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
source(file.path(dirname(if (length(here)) here[[1]] else "."), "_common.R"))

flags <- parse_flags()
module <- flag_str(flags, "module")
have <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# package -> role, the presets that rely on it, and whether it is core/required.
catalogue <- list(
  list(package = "GenomicRanges", role = "peak/interval container and overlap algebra", modules = 1:6, core = TRUE),
  list(package = "rtracklayer", role = "import/export BED, narrowPeak, broadPeak, GFF, bigWig", modules = 1:6, core = TRUE),
  list(package = "Rsamtools", role = "BAM access and index checks", modules = c(1, 2, 3, 6), core = TRUE),
  list(package = "GenomicAlignments", role = "summarizeOverlaps read counting and FRiP", modules = c(1, 3, 6), required = TRUE),
  list(package = "SummarizedExperiment", role = "assay() container for count matrices (transitive)", modules = c(1, 3, 6), core = TRUE),
  list(package = "S4Vectors", role = "mcols() metadata on regions (transitive)", modules = c(1, 2, 3), core = TRUE),
  list(package = "DiffBind", role = "occupancy/affinity analysis and dba.analyze", modules = c(1, 3), required = TRUE),
  list(package = "DESeq2", role = "negative-binomial differential binding (preset 1)", modules = c(1), required = TRUE),
  list(package = "edgeR", role = "QL differential binding (preset 3) and csaw engine", modules = c(2, 3), required = TRUE),
  list(package = "csaw", role = "window-based differential binding for broad marks (preset 2)", modules = c(2), required = TRUE),
  list(package = "ChIPseeker", role = "peak annotation, TSS profile and heatmap", modules = c(1, 3, 5), required = TRUE),
  list(package = "GenomicFeatures", role = "TxDb class infrastructure (transitive core); the GFF builder moved to txdbmaker and loadDb lives in AnnotationDbi", modules = c(1, 3, 5), required = TRUE),
  list(package = "txdbmaker", role = "makeTxDbFromGFF for user GTF/GFF3 annotation (split out of GenomicFeatures in Bioc >= 3.22)", modules = c(1, 3, 5), required = TRUE),
  list(package = "AnnotationDbi", role = "loadDb for a user prebuilt TxDb .sqlite (transitive)", modules = c(1, 3, 5), core = TRUE),
  list(package = "ComplexHeatmap", role = "peak-overlap UpSet plot", modules = c(4), required = TRUE),
  list(package = "ggplot2", role = "sample PCA and exploration figures", modules = c(6), required = TRUE),
  list(package = "ragg", role = "high-quality PNG figure device", modules = 1:6, required = FALSE),
  list(package = "ggupset", role = "ChIPseeker upsetplot fallback for a single annotated set (a Suggests, not shipped)", modules = c(1, 3, 5), required = FALSE),
  list(package = "TxDb.Hsapiens.UCSC.hg19.knownGene", role = "build-time ChIPseeker conda dependency; NOT a target-genome annotation source", modules = integer(0), required = FALSE),
  list(package = "idr", role = "formal IDR reproducibility test (pip)", modules = c(1, 3), required = FALSE),
  list(package = "profileplyr", role = "alternative coverage heatmap/profile engine", modules = c(1, 3, 5), required = FALSE),
  list(package = "org.Hs.eg.db", role = "human gene symbols for annotation", modules = c(1, 3, 5), required = FALSE),
  list(package = "org.Mm.eg.db", role = "mouse gene symbols for annotation", modules = c(1, 3, 5), required = FALSE)
)

probe <- function(entry) {
  installed <- have(entry$package)
  version <- if (installed) {
    tryCatch(as.character(utils::packageVersion(entry$package)), error = function(e) NA_character_)
  } else NA_character_
  list(
    package = entry$package,
    role = entry$role,
    installed = installed,
    version = version,
    core = isTRUE(entry$core),
    required = isTRUE(entry$required),
    modules = chr(as.integer(entry$modules))
  )
}

packages <- lapply(catalogue, probe)
absent <- function(entries) {
  vapply(Filter(function(entry) !have(entry$package), entries),
         function(entry) entry$package, character(1))
}

# A preset's core path needs only its required packages; optional ones stay
# unavailable here. core_ready avoids both overclaiming and refusing a preset
# whose main workflow works.
modules <- lapply(1:6, function(number) {
  needed <- Filter(function(entry) number %in% entry$modules, catalogue)
  required <- Filter(function(entry) isTRUE(entry$required), needed)
  optional <- Filter(function(entry) !isTRUE(entry$required), needed)
  missing_required <- absent(required)
  list(
    module = number,
    packages_checked = length(needed),
    missing_required = chr(missing_required),
    missing_optional = chr(absent(optional)),
    core_ready = length(missing_required) == 0
  )
})

selected <- if (!is.null(module)) {
  Filter(function(entry) entry$module %in% as.integer(module), modules)
} else modules

# Gaps stated as facts with a concrete next step. Package gaps are detected; the
# upstream and motif gaps are structural (this plugin is pure downstream).
gaps <- list()
add_gap <- function(capability, reason, next_step) {
  gaps[[length(gaps) + 1L]] <<- list(capability = capability, reason = reason, next_step = next_step)
}
if (!have("csaw")) {
  add_gap("broad-mark window differential binding (csaw)",
          "bioconda has no osx-arm64 build; csaw ships on linux-x64 / darwin-x64 only",
          "use the DiffBind route for broad marks, run on linux-x64 / darwin-x64, or request a project environment (env_request)")
}
if (!have("idr")) {
  add_gap("IDR reproducibility test",
          "the idr pip package is not shipped in this conda environment",
          "request a project environment (env_request); overlap fraction is a transparent proxy meanwhile")
}
if (!have("profileplyr")) {
  add_gap("profileplyr coverage heatmaps",
          "profileplyr is not part of the locked stack",
          "ChIPseeker plotAvgProf2/peakHeatmap covers the same view; request a project environment for profileplyr")
}
if (!have("org.Hs.eg.db") && !have("org.Mm.eg.db")) {
  add_gap("gene symbols in peak annotation",
          "lean annotation strategy bundles no org.db; annotation reports gene ids without symbols",
          "pass --org-db with a user-installed org.db, or request a project environment (env_request)")
}
add_gap("motif enrichment (HOMER / MEME)",
        "motif discovery is outside the R downstream stack",
        "run HOMER/MEME upstream or request a project environment (env_request)")
add_gap("upstream FASTQ -> BAM -> peak calling",
        "this plugin is pure downstream: it starts from called peaks / counts / BAM / bigWig",
        "ChIP: nf-core/chipseq wrapper; CUT&Tag: user SEACR run (no wrapper yet)")

installed_count <- sum(vapply(packages, function(entry) entry$installed, logical(1)))

emit(list(
  environment = paste0("R ", getRversion()),
  packages_total = length(packages),
  packages_installed = installed_count,
  packages = packages,
  modules = selected,
  core_ready_modules = chr(vapply(Filter(function(e) e$core_ready, selected),
                                  function(e) as.character(e$module), character(1))),
  blocked_modules = chr(vapply(Filter(function(e) !e$core_ready, selected),
                               function(e) as.character(e$module), character(1))),
  gaps = gaps,
  notes = chr(c(
    if (installed_count < length(packages))
      "the plugin environment is read-only and prebuilt: a missing package needs a project environment (env_request), not an in-environment install",
    "annotation is lean: the target-genome TxDb comes from a user GTF/GFF3 or TxDb .sqlite; the bundled hg19 TxDb is a ChIPseeker build dependency, never the default annotation source"
  ))
))