# Motif enrichment asks which TF binding motifs are over-represented in each
# cluster's marker regions versus background — the regulatory hypothesis layer of
# scATAC (which TFs may drive each cluster's open chromatin). It needs three things
# the environment does not bundle: the marker regions (from the markers step), a
# genome FASTA (genome_fasta; snapatac2's lazy download is ~845 MB from a slow host,
# so pass a local file), and a motif database (motif_db .meme, else snapatac2 cis_bp,
# itself a download). Missing any of them, the cell reports the gap and continues.
# This is motif enrichment, not chromVAR deviation scores (chromVAR is R/pip and lives
# in a project environment); a full TF-gene GRN is a snapatac2 extension via inject.
if not markers:
    analysis_notes.append("motif skipped: no marker regions; run the markers step first")
    info("motif skipped: no marker regions")
elif not genome_fasta:
    analysis_notes.append("motif skipped: no genome_fasta (pass a local FASTA; the snapatac2 download is ~845MB)")
    info("motif skipped: genome_fasta is null")
elif not os.path.exists(genome_fasta):
    analysis_notes.append("motif skipped: genome_fasta not found: " + genome_fasta)
    info("genome_fasta not found: " + genome_fasta)
else:
    try:
        if motif_db:
            if not os.path.exists(motif_db):
                raise ValueError("motif_db not found: " + motif_db)
            motifs = snap.read_motifs(motif_db)
        else:
            info("loading cis_bp motifs (snapatac2 cache download)")
            motifs = snap.datasets.cis_bp(unique=True)
        regions = {str(k): list(v) for k, v in markers.items() if v}
        if not regions:
            raise ValueError("no non-empty marker regions to test")
        enrichment = snap.tl.motif_enrichment(motifs, regions, genome_fasta, method=motif_method)
        os.makedirs(tables_dir, exist_ok=True)
        motif_results = {}
        for _cl, _df in enrichment.items():
            motif_results[str(_cl)] = int(_df.height)
            if _df.height:
                _df.write_csv(os.path.join(tables_dir, f"motif_{_cl}.csv"))
        analysis_steps.append("motif")
        print(f"motif enrichment ({motif_method}, {len(motifs)} motifs): {motif_results}")
    except Exception as error:
        analysis_notes.append("motif skipped: " + str(error))
        info("motif skipped: " + str(error))