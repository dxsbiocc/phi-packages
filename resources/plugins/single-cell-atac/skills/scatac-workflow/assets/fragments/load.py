# data_path, genome, chrom_sizes, min_num_fragments, is_paired and n_jobs come
# from the config cell. scATAC starts from FRAGMENTS (the aligned, deduplicated
# Tn5 insertion intervals), not from a gene count matrix: the features of record
# are genome bins (tiles) and, later, called peaks. import_fragments counts
# fragments per barcode and keeps the object file-backed so a large experiment
# never has to fit in RAM. chrom_sizes is bundled offline for hg38/mm10/hg19/mm39
# (resolve_chrom_sizes in setup); annotation and FASTA are not needed to import.
import shutil
import tempfile

kind = detect_kind(data_path)
loaded_kind = kind
info("input kind: " + kind)

# Work on a temporary backed file, never the user's input: import_fragments creates
# it, and a saved .h5ad is copied first so read-write backing cannot mutate the
# original. The save step exports the finished object to output_object.
_backing = tempfile.mktemp(suffix=".h5ad", prefix="scatac-")

if kind == "anndata":
    shutil.copyfile(data_path, _backing)
    adata = snap.read(_backing, backed="r+")
    genome_info = {"source": "saved-object", "input_kind": kind}
elif kind in ("10x-matrix-dir", "10x-h5"):
    adata = snap.read_10x_mtx(data_path)
    genome_info = {"source": "10x-matrix", "input_kind": kind}
    info("read a 10x peak/tile matrix; tiles and peaks are already the features")
elif kind in ("fragment", "fragment-like-table"):
    sizes, genome_info = resolve_chrom_sizes()
    adata = snap.pp.import_fragments(
        data_path, chrom_sizes=sizes, sorted_by_barcode=False,
        min_num_fragments=min_num_fragments, file=_backing,
        is_paired=is_paired, n_jobs=n_jobs,
    )
elif kind == "bam":
    # A BAM is turned into fragments first (Tn5 shift + barcode tag), then imported.
    frag = _backing + ".fragments.tsv.gz"
    snap.pp.make_fragment_file(data_path, frag, is_paired=is_paired, source="10x")
    sizes, genome_info = resolve_chrom_sizes()
    adata = snap.pp.import_fragments(
        frag, chrom_sizes=sizes, sorted_by_barcode=False,
        min_num_fragments=min_num_fragments, file=_backing,
        is_paired=is_paired, n_jobs=n_jobs,
    )
    genome_info["fragments_from_bam"] = frag
else:
    raise ValueError(
        "unsupported scATAC input (" + kind + "); expected a fragments.tsv.gz, a "
        "position-sorted BAM, a saved snapatac2 .h5ad, or a 10x peak/tile matrix "
        "directory/.h5. Raw FASTQ alignment belongs to a Cell Ranger ATAC wrapper."
    )

genome_info["input_kind"] = kind
print(f"loaded {adata.n_obs} barcodes x {adata.n_vars} features ({kind})")
# n_fragment is the per-cell fragment count snapatac2 imports; it is the primary
# QC metric (the scATAC analogue of total_counts). frac_dup / frac_mito come with it.
_nfrag = obs_col(adata, "n_fragment")
if _nfrag is not None and len(_nfrag):
    print(f"fragments per barcode: min {int(_nfrag.min())}, median {int(np.median(_nfrag))}, max {int(_nfrag.max())}")
if adata.n_obs == 0:
    raise ValueError(
        "no barcode survived import; lower min_num_fragments, check the genome build "
        "matches the fragments, or verify the file is a fragments table (chr/start/end/barcode/count)."
    )