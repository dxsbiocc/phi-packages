# Save the object so a later session, or a script, resumes from here instead of
# recomputing. The path is output_object in the config cell.
os.makedirs(os.path.dirname(os.path.abspath(output_object)), exist_ok=True)
adata.write_h5ad(output_object)
print(f"saved {adata.n_obs} cells x {adata.n_vars} genes -> {output_object}")