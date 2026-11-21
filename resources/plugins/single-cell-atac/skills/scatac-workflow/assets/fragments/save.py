# Save the finished object so a later session or script resumes from here instead of
# recomputing. The path is output_object in the config cell. snapatac2 uses .write()
# (not write_h5ad) to export a copy; the cluster labels, spectral/UMAP embeddings,
# QC metrics and any integrated embedding travel with the object. The peak and gene
# activity matrices are separate files written by their own steps (tables_dir).
os.makedirs(os.path.dirname(os.path.abspath(output_object)), exist_ok=True)
adata.write(output_object)
print(f"saved {adata.n_obs} cells x {adata.n_vars} features -> {output_object}")
print(f"steps run: {', '.join(analysis_steps) if analysis_steps else '(none recorded)'}")
if genome_info:
    print(f"genome: {genome_info.get('genome') or genome_info.get('source')}, input: {loaded_kind}")
for _obj in (peak_adata, gene_adata):
    if _obj is not None:
        try:
            _obj.close()
        except Exception:
            pass