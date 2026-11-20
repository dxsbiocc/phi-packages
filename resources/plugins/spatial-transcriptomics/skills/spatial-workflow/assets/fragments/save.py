# Save the object so a later session, or a script, resumes from here instead of
# recomputing. The path is output_object in the config cell. Spatial coordinates
# (obsm['spatial']), the raw counts layer and every spatial statistic written to
# uns/obsp travel with the object, so downstream steps can map results to tissue.
os.makedirs(os.path.dirname(os.path.abspath(output_object)), exist_ok=True)
adata.write_h5ad(output_object)
print(f"saved {adata.n_obs} spatial units x {adata.n_vars} genes -> {output_object}")
print(f"steps run: {', '.join(analysis_steps) if analysis_steps else '(none recorded)'}")
if platform_info:
    print(f"platform: {platform_info.get('platform')}, input: {platform_info.get('input_kind')}")