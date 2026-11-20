# Cell-level tests treat every cell as independent, which inflates significance
# because cells from one donor are not independent: the experimental unit is the
# donor, not the cell. Sum RAW counts per sample (x condition group) first, then
# test the samples with edgeR/DESeq2. Power comes from more donors, not more
# cells per donor.
import re

if sample_column not in adata.obs.columns:
    print(f"No '{sample_column}' column in this object; skipping pseudobulk.")
    print(f"Add one with adata.obs['{sample_column}'] = ... then re-run this cell.")
else:
    if group_column and group_column not in adata.obs.columns:
        raise ValueError(
            f"group column '{group_column}' not in obs; available: " + ", ".join(adata.obs.columns)
        )
    if "counts" not in adata.layers:
        raise ValueError("object has no counts layer; pseudobulk must start from raw counts")
    group_provided = bool(group_column)
    counts = adata.layers["counts"]
    sample_values = adata.obs[sample_column].astype(str).to_numpy()
    group_values = (
        adata.obs[group_column].astype(str).to_numpy()
        if group_provided
        else np.full(adata.n_obs, "all")
    )
    unit = (
        np.array([f"{s}|{g}" for s, g in zip(sample_values, group_values)])
        if group_provided
        else sample_values
    )
    unit_levels = list(dict.fromkeys(unit))
    per_unit_cells = pd.Series(unit).value_counts()
    small_units = [name for name, count in per_unit_cells.items() if count < min_cells]

    columns = {}
    design = []
    for name in unit_levels:
        mask = unit == name
        cells = int(mask.sum())
        if cells < min_cells:
            info(f"skipping {name}: {cells} cells < {min_cells}")
            continue
        summed = np.asarray(counts[mask].sum(axis=0)).ravel()
        parts = name.split("|")
        identifier = re.sub(r"[^A-Za-z0-9._-]", "_", name)
        columns[identifier] = summed
        design.append(
            {
                "column": identifier,
                "sample": parts[0],
                "group": parts[1] if group_provided else "all",
                "cells": cells,
            }
        )
    if not columns:
        raise ValueError(
            f"no unit has at least {min_cells} cells; lower --min-cells or check "
            f"--sample-column ({sample_column})"
        )

    matrix_out = pd.DataFrame(columns, index=adata.var_names)
    # The filter counts pseudobulk columns (samples), not cells: a gene must be
    # observed in enough independent samples to be testable at all.
    expressed = (matrix_out >= min_count).sum(axis=1) >= min_cells_expressing
    dropped = int((~expressed).sum())
    matrix_out = matrix_out[expressed]
    if matrix_out.shape[0] == 0:
        raise ValueError(
            f"no gene passed the expression filter (--min-count {min_count}, "
            f"--min-cells-expressing {min_cells_expressing})"
        )

    condition_column = group_column if group_provided else "group"
    sample_sheet = pd.DataFrame(
        {
            "sample": [entry["column"] for entry in design],
            "source_sample": [entry["sample"] for entry in design],
            condition_column: [entry["group"] for entry in design],
            "cells": [entry["cells"] for entry in design],
        }
    )

    written = []
    if counts_output:
        os.makedirs(os.path.dirname(os.path.abspath(counts_output)), exist_ok=True)
        matrix_out.reset_index(names="gene_id").to_csv(counts_output, sep="\t", index=False)
        info("counts " + counts_output)
        written.append(counts_output)
    if samples_output:
        os.makedirs(os.path.dirname(os.path.abspath(samples_output)), exist_ok=True)
        sample_sheet.to_csv(samples_output, sep="\t", index=False)
        info("sample sheet " + samples_output)
        written.append(samples_output)

    group_counts = sample_sheet[condition_column].value_counts()
    notes = []
    if group_provided:
        thin = [name for name, count in group_counts.items() if count < min_samples_per_group]
        if thin:
            notes.append(
                f"group(s) with fewer than {min_samples_per_group} samples: "
                + ", ".join(thin)
                + "; a two-group test needs biological replicates"
            )
        if len(group_counts) < 2:
            notes.append("only one group is present, so no comparison can be made")
    else:
        notes.append("no --group-column: columns are one per sample and carry no condition")
    if small_units:
        notes.append(f"skipped {len(small_units)} unit(s) under {min_cells} cells")
    if dropped:
        notes.append(f"dropped {dropped} gene(s) below the expression filter")

    pb_info = {
        "outputs": written,
        "units": int(matrix_out.shape[1]),
        "genes": int(matrix_out.shape[0]),
        "design": design,
        "sample_column": sample_column,
        "group_column": group_column if group_provided else None,
        "wrapper_hint": {
            "counts": counts_output,
            "samplesheet": samples_output,
            "sample_id_col": "sample",
            "gene_id_col": "gene_id",
            "contrast_variable": condition_column,
        }
        if group_provided
        else None,
        "parameters": {
            "min_cells": min_cells,
            "min_count": min_count,
            "min_cells_expressing": min_cells_expressing,
        },
        "dropped_genes": dropped,
        "notes": notes,
    }
    print(f"pseudobulk: {pb_info['genes']} genes x {pb_info['units']} units")
    if written:
        print("Next: run the deseq2 wrapper with counts=" + str(counts_output)
              + " and samplesheet=" + str(samples_output))