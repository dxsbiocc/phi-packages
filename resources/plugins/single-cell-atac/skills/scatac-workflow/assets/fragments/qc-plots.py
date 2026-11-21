# Look at the distributions BEFORE choosing thresholds. A threshold is a biological
# decision read off these plots, not a round number accepted blindly. plot_tag is
# "before" on the first pass and "after" once filtering has run, so the two panels
# can be compared side by side. Every figure is written as interactive HTML.
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# The nucleosome ladder: a clean ATAC library peaks below ~100 bp (nucleosome-free)
# with ~200/400/600 bp shoulders. A flat ladder means a poor library.
if "frag_size_distr" in qc_scored:
    try:
        show_plot(snap.pl.frag_size_distr(adata, show=False), f"qc-frag-size-{plot_tag}")
    except Exception as error:
        info("frag_size_distr plot skipped: " + str(error))

# TSS enrichment vs fragment count: high-TSSe cells are the good libraries.
if "tsse" in qc_scored:
    try:
        show_plot(snap.pl.tsse(adata, show=False), f"qc-tsse-{plot_tag}")
    except Exception as error:
        info("tsse plot skipped: " + str(error))

# Per-cell import metrics. n_fragment spans orders of magnitude, so a log x-axis
# reveals the empty-droplet shoulder that a linear axis hides.
_metrics = [m for m in ("n_fragment", "frac_dup", "frac_mito") if m in qc_scored]
if _metrics:
    fig = make_subplots(rows=1, cols=len(_metrics), subplot_titles=_metrics)
    for _i, _m in enumerate(_metrics, start=1):
        _v = obs_col(adata, _m)
        fig.add_trace(go.Histogram(x=_v, nbinsx=60, name=_m), row=1, col=_i)
        if _m == "n_fragment":
            fig.update_xaxes(type="log", row=1, col=_i)
    fig.update_layout(
        title_text=f"scATAC QC distributions ({plot_tag})",
        showlegend=False, height=380, bargap=0.02,
    )
    show_plot(fig, f"qc-metrics-{plot_tag}")
print(f"QC plots written ({plot_tag}): {', '.join(k for k in figure_paths if k.endswith(plot_tag))}")