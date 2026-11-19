# Doublets are two cells in one droplet: they look like a novel "hybrid" state
# and can invent whole clusters. Scrublet simulates doublets from the observed
# counts and scores each cell against them. Run it per sample when samples were
# loaded separately; the expected rate scales with how many cells were loaded.
if not detect_doublets:
    print("doublet detection skipped (detect_doublets is off)")
elif not optional_module("scrublet"):
    doublet_result = {
        "method": "none",
        "error": "scrublet is not installed in this environment; "
        "install it into a project environment to enable doublet detection",
    }
    print(doublet_result["error"])
else:
    try:
        rate = doublet_expected_rate
        if rate is None:
            rate = max(0.01, min(0.15, 0.008 * (adata.n_obs / 1000)))
        sc.pp.scrublet(adata, expected_doublet_rate=rate, rng=123, use="layers.counts")
        predicted = adata.obs["predicted_doublet"].astype(bool)
        doublet_result = {
            "method": "scrublet",
            "expected_rate": round(rate, 4),
            "singlet": int((~predicted).sum()),
            "doublet": int(predicted.sum()),
        }
        print(f"scrublet flagged {doublet_result['doublet']} predicted doublets "
              f"({doublet_result['doublet'] / adata.n_obs:.1%})")
    except Exception as error:  # scrublet can fail on tiny or sparse inputs
        info("scrublet failed: " + str(error))
        doublet_result = {"method": "scrublet", "error": str(error)}
# The filter cell below removes the flagged cells; scoring alone changes nothing.