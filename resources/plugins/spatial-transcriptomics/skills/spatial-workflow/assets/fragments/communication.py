# Cell-cell communication infers which ligand-receptor pairs could link two
# domains: a ligand expressed in domain A and its receptor in domain B. squidpy's
# ligrec runs a CellPhoneDB-style permutation test — it shuffles cluster labels
# n_perms times and asks how often the mean ligand/receptor expression exceeds the
# observed value, giving a p-value per pair per cluster combination.
#
# ligrec normally fetches interactions from OmniPath online. To run offline in a
# locked environment we pass a curated ligand-receptor DataFrame (source=ligand,
# target=receptor). The built-in set below is a small starter covering major
# signalling families — NOT a comprehensive database. For full coverage pass
# interactions_path (a TSV/CSV with source/target or ligand/receptor columns), or
# use liana/CellPhoneDB in a project environment. Gene symbols are matched
# case-insensitively (ligrec uppercases), so a human-symbol resource also works on
# mouse data by orthologous symbol.
import squidpy as sq
import matplotlib.pyplot as plt

if group_by not in adata.obs.columns:
    raise ValueError(f"group column '{group_by}' not in obs")
if not isinstance(adata.obs[group_by].dtype, pd.CategoricalDtype):
    adata.obs[group_by] = adata.obs[group_by].astype("category")
if len(adata.obs[group_by].cat.categories) < 2:
    info("communication needs >= 2 groups; skipping")
    analysis_notes.append("communication skipped: fewer than 2 groups")
    ccc_results = {"skipped": "fewer than 2 groups"}
else:
    if interactions_path:
        if not os.path.exists(interactions_path):
            raise ValueError("interactions file not found: " + interactions_path)
        sep = "," if interactions_path.lower().endswith(".csv") else "\t"
        lr = pd.read_csv(interactions_path, sep=sep)
        if {"source", "target"} <= set(lr.columns):
            lr = lr[["source", "target"]]
        elif {"ligand", "receptor"} <= set(lr.columns):
            lr = lr.rename(columns={"ligand": "source", "receptor": "target"})[["source", "target"]]
        else:
            raise ValueError("interactions file needs source/target or ligand/receptor columns")
        resource = os.path.basename(interactions_path)
    else:
        pairs = [
            # growth factors
            ("TGFB1", "TGFBR1"), ("TGFB1", "TGFBR2"), ("TGFB2", "TGFBR1"), ("BMP4", "BMPR1A"),
            ("BMP4", "BMPR2"), ("FGF2", "FGFR1"), ("FGF2", "FGFR2"), ("VEGFA", "FLT1"),
            ("VEGFA", "KDR"), ("PDGFB", "PDGFRB"), ("EGF", "EGFR"), ("HGF", "MET"),
            ("IGF1", "IGF1R"), ("NGF", "NGFR"), ("BDNF", "NTRK2"), ("INHBA", "ACVR1B"),
            # cytokines
            ("TNF", "TNFRSF1A"), ("TNF", "TNFRSF1B"), ("IL6", "IL6R"), ("IL6", "IL6ST"),
            ("IL1B", "IL1R1"), ("IL10", "IL10RA"), ("IFNG", "IFNGR1"), ("IFNG", "IFNGR2"),
            ("IL4", "IL4R"), ("IL13", "IL13RA1"), ("IL17A", "IL17RA"), ("IL2", "IL2RA"),
            ("IL15", "IL15RA"), ("IL18", "IL18R1"), ("IL33", "IL1RL1"),
            # chemokines
            ("CCL2", "CCR2"), ("CCL5", "CCR5"), ("CXCL12", "CXCR4"), ("CXCL8", "CXCR1"),
            ("CXCL8", "CXCR2"), ("CXCL10", "CXCR3"), ("CCL19", "CCR7"), ("CCL21", "CCR7"),
            ("CX3CL1", "CX3CR1"), ("CXCL13", "CXCR5"),
            # ECM / adhesion
            ("COL1A1", "ITGB1"), ("COL1A1", "ITGA2"), ("FN1", "ITGB1"), ("FN1", "ITGAV"),
            ("LAMA1", "ITGB1"), ("VCAM1", "ITGA4"), ("ICAM1", "ITGAL"), ("SPP1", "CD44"),
            ("SPP1", "ITGAV"), ("CDH1", "CDH1"), ("NCAM1", "NCAM1"),
            # Notch / Wnt
            ("DLL1", "NOTCH1"), ("DLL4", "NOTCH1"), ("JAG1", "NOTCH1"), ("JAG2", "NOTCH2"),
            ("WNT3A", "FZD1"), ("WNT5A", "FZD5"), ("WNT7A", "FZD9"), ("RSPO3", "LGR5"),
            # immune checkpoints / TNF superfamily
            ("CD40LG", "CD40"), ("FASLG", "FAS"), ("CD70", "CD27"), ("CD274", "PDCD1"),
            ("PDCD1LG2", "PDCD1"), ("CTLA4", "CD80"), ("CTLA4", "CD86"), ("CD80", "CD28"),
            ("CD86", "CD28"), ("LGALS9", "HAVCR2"), ("TNFSF10", "TNFRSF10A"), ("TNFSF11", "TNFRSF11A"),
            # complement / vascular / metabolic
            ("C3", "C3AR1"), ("C5", "C5AR1"), ("ANGPT1", "TEK"), ("ANGPT2", "TEK"),
            ("EDN1", "EDNRA"), ("APLN", "APLNR"), ("LEP", "LEPR"), ("ADIPOQ", "ADIPOR1"),
            ("GH1", "GHR"), ("PRL", "PRLR"),
        ]
        lr = pd.DataFrame(pairs, columns=["source", "target"])
        resource = "curated starter (inline)"

    # ligrec reads adata.X (no layer argument) and defaults to use_raw=True. Build a
    # minimal object whose X is the un-scaled lognorm snapshot, so the test runs on
    # log expression rather than z-scores, without copying obsm/layers.
    expr = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    ccc = sc.AnnData(X=expr, obs=adata.obs[[group_by]].copy(), var=adata.var.copy())
    # A screen, not a guarantee: if none of the resource's ligand/receptor genes are
    # present (wrong species/identifier type, or a panel that lacks them), ligrec has
    # nothing to test. Report the gap and continue rather than crash the notebook.
    resource_genes = set(lr["source"]) | set(lr["target"])
    if not (resource_genes & set(ccc.var_names)):
        info("communication skipped: no ligand-receptor resource genes present in this object")
        ccc_results = {"skipped": "no ligand-receptor resource genes present", "resource": resource}
        analysis_notes.append(
            "communication skipped: no resource genes present (check species/identifier type or pass interactions_path)"
        )
        result = None
    else:
        try:
            result = sq.gr.ligrec(
                ccc,
                cluster_key=group_by,
                interactions=lr,
                use_raw=False,
                n_perms=n_perms,
                threshold=0.01,
                corr_method="fdr_bh",
                corr_axis="clusters",
                n_jobs=1,
                show_progress_bar=False,
                seed=seed,
                copy=True,
            )
        except ValueError as error:
            if "no interactions remain" not in str(error):
                raise
            info(f"communication skipped: {error}")
            ccc_results = {"skipped": "no interactions remain after filtering by genes", "resource": resource}
            analysis_notes.append(
                "communication skipped: resource genes present but none passed ligrec's expression filter"
            )
            result = None

    if result is not None:
        adata.uns[f"{group_by}_ligrec"] = result
        analysis_steps.append("communication")

        means = result["means"]
        pvalues = result["pvalues"]
        # Stack every column level so the result is a Series indexed by
        # (ligand, receptor, source_cluster, target_cluster) across pandas versions.
        col_levels = list(range(means.columns.nlevels))
        means_long = means.stack(level=col_levels)
        pvals_long = pvalues.stack(level=col_levels)
        significant = means_long[(pvals_long < ccc_threshold) & (means_long > 0)]
        top = significant.sort_values(ascending=False).head(10)
        ccc_results = {
            "method": "ligrec",
            "resource": resource,
            "n_interactions_tested": int(len(means)),
            "n_significant": int(len(significant)),
            "pvalue_threshold": ccc_threshold,
            "n_perms": n_perms,
            "top": [
                {
                    "ligand": idx[0],
                    "receptor": idx[1],
                    "source_cluster": idx[2],
                    "target_cluster": idx[3],
                    "mean": round(float(val), 4),
                }
                for idx, val in top.items()
            ],
        }
        print(f"ligrec: {len(significant)} significant interactions (p < {ccc_threshold}) of {len(means)} tested")
        for entry in ccc_results["top"][:6]:
            print(f"  {entry['ligand']}->{entry['receptor']}  {entry['source_cluster']}|{entry['target_cluster']}  mean {entry['mean']}")

        try:
            sq.pl.ligrec(
                adata,
                cluster_key=group_by,
                pvalue_threshold=ccc_threshold,
                alpha=ccc_threshold,
                swap_axes=True,
            )
            show_plot(plt.gcf(), "communication-ligrec")
        except ValueError as error:
            info(f"ligrec plot skipped: {error}")

        # The raw ligrec dict holds MultiIndex-column DataFrames that h5ad cannot
        # serialize (tuple keys). After plotting, replace it with a tidy long table
        # with string columns so the saved object stays writable.
        tidy = means_long.rename("mean").to_frame()
        tidy["pvalue"] = pvals_long
        tidy["mean"] = np.asarray(tidy["mean"], dtype=float)
        tidy["pvalue"] = np.asarray(tidy["pvalue"], dtype=float)
        tidy.index.names = ["ligand", "receptor", "source_cluster", "target_cluster"]
        adata.uns[f"{group_by}_ligrec"] = tidy.reset_index()