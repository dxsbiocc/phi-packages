# Run the notebook fragment library as the single source of analysis code.
# The script tools are thin adapters: they translate CLI flags into the same
# config variables the notebook config cell defines, run the same fragments in
# order, and let scripts/reports/*.R serialize the bookkeeping into JSON.
# Nothing analytical lives here or in the adapters.

fragment_path <- function(step) {
  file.path(script_dir(), "..", "assets", "fragments", paste0(step, ".R"))
}

run_fragments <- function(steps) {
  for (step in steps) {
    path <- fragment_path(step)
    if (!file.exists(path)) fail("unknown fragment: ", step)
    info("fragment ", step)
    sys.source(path, envir = globalenv())
  }
  invisible(steps)
}