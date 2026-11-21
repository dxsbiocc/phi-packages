# Map CLI flags to the notebook config variables. This mirrors the config cell
# the renderer builds from assets/steps.json, so a script tool and a notebook run
# the same fragments with the same values. Sourced after _common.R (which provides
# parse_flags and the flag_* helpers) and before run_fragments().
flags <- parse_flags()

data_path <- flag_req(flags, "input")
project_name <- flag_str(flags, "project", "chip")
output_object <- flag_str(flags, "output")
figures_dir <- flag_str(flags, "figures_dir", "figures")
seed <- flag_int(flags, "seed", 42L)
sample_column <- flag_str(flags, "sample_column", "SampleID")

peak_paths <- flag_vec(flags, "peaks")
count_path <- flag_str(flags, "counts")
bam_paths <- flag_vec(flags, "bam")
bam_dir <- flag_str(flags, "bam_dir")
sample_sheet_path <- flag_str(flags, "sample_sheet")
paired_end <- flag_bool(flags, "paired_end", FALSE)
blacklist_path <- flag_str(flags, "blacklist")
genome <- flag_str(flags, "genome")

min_occupancy <- flag_int(flags, "min_occupancy", 2L)
normalize_method <- flag_str(flags, "normalize_method", "libsize")
contrast <- flag_vec(flags, "contrast")
min_members <- flag_int(flags, "min_members", 2L)
fdr_threshold <- flag_num(flags, "fdr", 0.05)
lfc_threshold <- flag_num(flags, "lfc", 0)
top_n <- flag_int(flags, "top_n", 20L)

csaw_width <- flag_int(flags, "csaw_width", 150L)
csaw_spacing <- flag_int(flags, "csaw_spacing", 50L)
csaw_filter <- flag_str(flags, "csaw_filter", "global")
fragment_size <- flag_int(flags, "fragment_size", 200L)

gtf_path <- flag_str(flags, "gtf")
txdb_path <- flag_str(flags, "txdb")
org_db <- flag_str(flags, "org_db")
tss_upstream <- flag_int(flags, "tss_upstream", 3000L)
tss_downstream <- flag_int(flags, "tss_downstream", 3000L)

# If the primary input is itself a DiffBind sample sheet, treat it as such so the
# load-peaks / load-bam fragments derive their paths from its columns.
if ((is.null(sample_sheet_path) || !nzchar(sample_sheet_path)) &&
    grepl("[.]csv$", data_path, ignore.case = TRUE) && file.exists(data_path)) {
  header <- tolower(readLines(data_path, n = 1))
  if (grepl("sampleid", header) && grepl("peaks", header)) sample_sheet_path <- data_path
}