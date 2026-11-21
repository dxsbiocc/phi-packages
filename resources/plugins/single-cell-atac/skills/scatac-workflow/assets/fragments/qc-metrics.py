# scATAC QC rests on three per-cell metrics that import_fragments already wrote:
#   n_fragment   fragments per barcode    (low = empty droplet / failed library)
#   frac_dup     duplication fraction     (high = over-sequenced, low complexity)
#   frac_mito    mitochondrial fraction   (high = damaged / dying cell)
# Two more need a computation:
#   frag_size_distr  the fragment-length histogram — a clean ATAC library shows the
#                    nucleosome ladder (~<100 bp nucleosome-free, ~200/400/600 bp
#                    mono-/di-/tri-nucleosome). A flat or single-peak ladder is a
#                    library-quality red flag, not a threshold to tune.
#   tsse         TSS enrichment — signal at transcription start sites relative to
#                flanks; the single best signal-to-noise metric for ATAC. It needs a
#                gene annotation, so it is skipped (with a gap message) when gene_anno
#                is null rather than triggering a slow ~50 MB download.
qc_scored = []
for _m in ("n_fragment", "frac_dup", "frac_mito"):
    _v = obs_col(adata, _m)
    if _v is not None and len(_v):
        qc_scored.append(_m)
        print(f"{_m}: min {np.min(_v):.4g}, median {np.median(_v):.4g}, max {np.max(_v):.4g}")

# Fragment size distribution is offline (it reads the imported fragments).
try:
    snap.metrics.frag_size_distr(adata, inplace=True)
    qc_scored.append("frag_size_distr")
except Exception as error:
    analysis_notes.append("frag_size_distr unavailable: " + str(error))
    info("frag_size_distr skipped: " + str(error))

# TSS enrichment needs a gene annotation. Without one, report the gap and continue.
if gene_anno:
    if not os.path.exists(gene_anno):
        analysis_notes.append("gene_anno path not found: " + gene_anno + "; tsse skipped")
        info("gene_anno not found: " + gene_anno)
    else:
        try:
            snap.metrics.tsse(adata, gene_anno, inplace=True)
            qc_scored.append("tsse")
            _t = obs_col(adata, "tsse")
            if _t is not None and len(_t):
                print(f"tsse: min {np.min(_t):.3f}, median {np.median(_t):.3f}, max {np.max(_t):.3f}")
        except Exception as error:
            analysis_notes.append("tsse failed: " + str(error))
            info("tsse failed: " + str(error))
else:
    analysis_notes.append("tsse skipped: no gene_anno (pass a GTF/GFF3 to enable TSS enrichment)")
    info("tsse skipped: gene_anno is null")

# Snapshot the object BEFORE any filtering. snapatac2 filters in place, so a
# reference to adata would shrink too; capture plain numbers and per-sample counts.
def _per_sample(obj):
    col = obs_col(obj, sample_column)
    if col is None:
        return None
    names, counts = np.unique(col.astype(str), return_counts=True)
    return [{"sample": str(n), "cells": int(c)} for n, c in zip(names, counts)]


obj_before = {"cells": adata.n_obs, "genes": adata.n_vars, "per_sample": _per_sample(adata)}
print(f"QC metrics computed: {', '.join(qc_scored)}")