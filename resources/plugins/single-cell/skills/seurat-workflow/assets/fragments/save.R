# Save the object so a later session, or a script, resumes from here instead of
# recomputing. The path is output_object in the config cell.
dir.create(dirname(output_object), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, output_object)
cat("saved", output_object, "\n")