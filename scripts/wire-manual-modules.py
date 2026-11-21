#!/usr/bin/env python3
"""One-off manual wiring for the 27 modules the generator skipped (exotic stageAs)."""
import json, os, re, sys
import yaml

SRC = '/Users/dengxsh/Downloads/Work/modules-master/modules/nf-core'
DST = 'resources/wrappers/modules/nf-core'
TD = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/'

RESOURCES = {
    'process_single': (1, '2.GB', '1.h'),
    'process_low': (2, '4.GB', '2.h'),
    'process_medium': (4, '8.GB', '4.h'),
    'process_high': (8, '16.GB', '8.h'),
    'process_long': (2, '4.GB', '20.h'),
    'process_high_memory': (8, '64.GB', '8.h'),
}

HEADER = """#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { %(proc)s } from '../main.nf'

"""

CONFIG = """// Launcher config for the %(proc)s wrapper adapter. Not part of the fixed
// wrapper/{main.nf,params.json,wrapper.yaml} contract — Nextflow picks this
// up automatically because it sits beside wrapper/main.nf.

process {
    cpus   = %(cpus)d
    memory = '%(mem)s'
    time   = '%(time)s'

    withName: '%(proc)s' {
        publishDir = [
            path: { "${params.outdir}/%(outdir_name)s" },
            mode: 'copy'
        ]
%(extargs)s    }
}

profiles {
    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
        conda.enabled = false
        singularity.enabled = false
    }
    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
        docker.enabled = false
        conda.enabled = false
    }
    conda {
        conda.enabled = true
        docker.enabled = false
        singularity.enabled = false
    }
}
"""

SPECS = []

def spec(rel, proc, name, summary, params, defaults, body, label, args=True, note=None):
    SPECS.append(dict(rel=rel, proc=proc, name=name, summary=summary, params=params,
                      defaults=defaults, body=body, label=label, args=args, note=note))

def P(kind, typ, desc, required=False, enum=None):
    d = {'kind': kind, 'type': typ, 'required': required, 'description': desc}
    if enum: d['enum'] = enum
    return d

OUTDIR = P('output', 'path', 'Output directory for the results.', True)

# ---------------------------------------------------------------- ashlar
spec('ashlar', 'ASHLAR', 'ashlar', 'Stitch and register multi-cycle microscopy images with ASHLAR.',
     {'images': P('input', 'path', 'Directory of per-cycle OME-TIFF images.', True),
      'opt_dfp': P('input', 'file', 'Optional dark-field profile image (dfp).'),
      'opt_ffp': P('input', 'file', 'Optional flat-field profile image (ffp).'),
      'outdir': OUTDIR},
     {'images': None,
      'opt_dfp': TD + 'imaging/ome-tiff/cycif-tonsil-dfp.ome.tif',
      'opt_ffp': TD + 'imaging/ome-tiff/cycif-tonsil-ffp.ome.tif',
      'outdir': 'results'},
     """workflow {
    images_ch = Channel.fromPath("${params.images}/*").collect()
    def dfp = params.opt_dfp ? file(params.opt_dfp, checkIfExists: true) : []
    def ffp = params.opt_ffp ? file(params.opt_ffp, checkIfExists: true) : []

    ASHLAR(images_ch.map { imgs -> [[id: 'ashlar'], imgs] }, dfp, ffp)
}
""", 'process_single')

# ---------------------------------------------------------------- bases2fastq
spec('bases2fastq', 'BASES2FASTQ', 'bases2fastq', 'Demultiplex Element AVITI/Polony Bases runs into FASTQ with bases2fastq.',
     {'run_manifest': P('input', 'file', 'RunManifest.csv describing the run.', True),
      'run_dir': P('input', 'path', 'Run directory produced by the instrument (staged as input_dir).', True),
      'outdir': OUTDIR},
     {'run_manifest': 'https://raw.githubusercontent.com/nf-core/test-datasets/demultiplex/testdata/sim-data/RunManifest.csv',
      'run_dir': None, 'outdir': 'results'},
     """workflow {
    BASES2FASTQ(Channel.value([
        [id: 'run'],
        file(params.run_manifest, checkIfExists: true),
        file(params.run_dir, checkIfExists: true)
    ]))
}
""", 'process_high')

# ---------------------------------------------------------------- cellranger/count
spec('cellranger/count', 'CELLRANGER_COUNT', 'cellranger count', 'Run 10x Genomics Cell Ranger count on GEX FASTQs against a Cell Ranger reference.',
     {'fastqs': P('input', 'path', 'Directory of GEX FASTQ files (*.fastq.gz).', True),
      'reference': P('input', 'path', 'Cell Ranger reference directory (from cellranger mkref).', True),
      'outdir': OUTDIR},
     {'fastqs': None, 'reference': None, 'outdir': 'results'},
     """workflow {
    reads_ch = Channel.fromPath("${params.fastqs}/*.fastq.gz").collect()
    ref_ch   = Channel.value([[id: 'reference'], file(params.reference, checkIfExists: true)])

    CELLRANGER_COUNT(reads_ch.map { r -> [[id: 'cellranger', single_end: false], r] }, ref_ch)
}
""", 'process_high')

# ---------------------------------------------------------------- cellranger/multi
spec('cellranger/multi', 'CELLRANGER_MULTI', 'cellranger multi', 'Run 10x Genomics Cell Ranger multi (GEX plus optional VDJ/antibody/BEAM/CMO/CRISPR libraries).',
     {'gex_fastqs': P('input', 'path', 'Directory of gene-expression FASTQs.', True),
      'gex_reference': P('input', 'path', 'GEX reference directory (Cell Ranger transcriptome).', True),
      'vdj_fastqs': P('input', 'path', 'Optional VDJ FASTQ directory.'),
      'ab_fastqs': P('input', 'path', 'Optional antibody-capture FASTQ directory.'),
      'beam_fastqs': P('input', 'path', 'Optional BEAM FASTQ directory.'),
      'cmo_fastqs': P('input', 'path', 'Optional CMO FASTQ directory.'),
      'crispr_fastqs': P('input', 'path', 'Optional CRISPR FASTQ directory.'),
      'gex_frna_probeset': P('input', 'file', 'Optional fixed-RNA probe set CSV.'),
      'gex_targetpanel': P('input', 'file', 'Optional target-panel CSV.'),
      'vdj_reference': P('input', 'path', 'Optional VDJ reference directory.'),
      'vdj_primer_index': P('input', 'file', 'Optional VDJ primer index FASTA.'),
      'fb_reference': P('input', 'path', 'Optional feature-barcode reference CSV.'),
      'beam_antigen_panel': P('input', 'file', 'Optional BEAM antigen panel CSV.'),
      'beam_control_panel': P('input', 'file', 'Optional BEAM control panel CSV.'),
      'cmo_reference': P('input', 'path', 'Optional CMO reference CSV.'),
      'cmo_barcodes': P('input', 'file', 'Optional CMO barcodes CSV.'),
      'cmo_barcode_assignment': P('input', 'file', 'Optional CMO sample barcode assignment CSV.'),
      'frna_sampleinfo': P('input', 'file', 'Optional fixed-RNA sample info CSV.'),
      'ocm_barcodes': P('input', 'file', 'Optional OCM barcodes CSV.'),
      'skip_renaming': P('option', 'boolean', 'Skip Cell Ranger sample renaming.'),
      'outdir': OUTDIR},
     {'gex_fastqs': None, 'gex_reference': None, 'vdj_fastqs': None, 'ab_fastqs': None,
      'beam_fastqs': None, 'cmo_fastqs': None, 'crispr_fastqs': None,
      'gex_frna_probeset': None, 'gex_targetpanel': None, 'vdj_reference': None,
      'vdj_primer_index': None, 'fb_reference': None, 'beam_antigen_panel': None,
      'beam_control_panel': None, 'cmo_reference': None, 'cmo_barcodes': None,
      'cmo_barcode_assignment': None, 'frna_sampleinfo': None, 'ocm_barcodes': None,
      'skip_renaming': False, 'outdir': 'results'},
     """// REVIEW: only the GEX-only and full-library combinations are exercised upstream;
// other partial combinations are not runtime-verified.
workflow {
    def fq  = { p -> p ? file(p, checkIfExists: true).listFiles().findAll { it.isFile() }.toList() : [] }
    def ref = { p -> p ? file(p, checkIfExists: true) : [] }

    CELLRANGER_MULTI(
        [id: 'multi'],
        Channel.value([[id: 'gex'],    fq(params.gex_fastqs),    '']),
        Channel.value([[id: 'vdj'],    fq(params.vdj_fastqs),    '']),
        Channel.value([[id: 'ab'],     fq(params.ab_fastqs),     '']),
        Channel.value([[id: 'beam'],   fq(params.beam_fastqs),   '']),
        Channel.value([[id: 'cmo'],    fq(params.cmo_fastqs),    '']),
        Channel.value([[id: 'crispr'], fq(params.crispr_fastqs), '']),
        Channel.value([[id: 'gex_reference'], ref(params.gex_reference)]),
        ref(params.gex_frna_probeset),
        ref(params.gex_targetpanel),
        ref(params.vdj_reference),
        ref(params.vdj_primer_index),
        ref(params.fb_reference),
        ref(params.beam_antigen_panel),
        ref(params.beam_control_panel),
        ref(params.cmo_reference),
        ref(params.cmo_barcodes),
        ref(params.cmo_barcode_assignment),
        ref(params.frna_sampleinfo),
        ref(params.ocm_barcodes),
        params.skip_renaming
    )
}
""", 'process_high')

# ---------------------------------------------------------------- comebin/runcomebin
spec('comebin/runcomebin', 'COMEBIN_RUNCOMEBIN', 'comebin runcomebin', 'Bin metagenomic contigs with COMEBIN using co-abundance across samples.',
     {'assembly': P('input', 'file', 'Assembly contigs FASTA (optionally gzipped).', True),
      'bams': P('input', 'path', 'Directory of per-sample sorted BAMs (*.bam).', True),
      'outdir': OUTDIR},
     {'assembly': TD + 'genomics/prokaryotes/metagenome/binning/contigs.fasta.gz',
      'bams': None, 'outdir': 'results'},
     """workflow {
    bams_ch = Channel.fromPath("${params.bams}/*.bam").collect()

    COMEBIN_RUNCOMEBIN(bams_ch.map { b ->
        [[id: 'comebin'], file(params.assembly, checkIfExists: true), b]
    })
}
""", 'process_high')

# ---------------------------------------------------------------- diann
spec('diann', 'DIANN', 'DIA-NN', 'Quantify DIA proteomics runs with DIA-NN (raw files, spectral-library prediction, or --use-quant).',
     {'ms_files': P('input', 'path', 'Directory of raw MS files (.raw/.mzML/.dia). Omit for library prediction or quant reuse.'),
      'ms_file_names': P('option', 'string', 'Comma-separated MS file basenames to use with quant (instead of staging raw files).'),
      'fasta': P('input', 'file', 'FASTA proteome for spectral-library prediction.'),
      'library': P('input', 'file', 'Existing spectral library (.speclib).'),
      'quant': P('input', 'path', 'Directory of preprocessed .quant files (enables --use-quant).'),
      'outdir': OUTDIR},
     {'ms_files': None, 'ms_file_names': None,
      'fasta': TD + 'proteomics/diann/REF_EColi_K12_UPS1_combined_subset_100.fasta',
      'library': None, 'quant': None, 'outdir': 'results'},
     """// The vendored module refuses -profile conda; use docker/singularity.
workflow {
    def ms = params.ms_files
        ? file(params.ms_files, checkIfExists: true).listFiles().findAll { it.isFile() }.toList()
        : []
    def names = params.ms_file_names ? params.ms_file_names.tokenize(',')*.trim() : []
    def quant = params.quant
        ? file(params.quant, checkIfExists: true).listFiles().findAll { it.isFile() }.toList()
        : []

    DIANN(Channel.value([
        [id: 'diann'], ms, names,
        params.fasta   ? file(params.fasta, checkIfExists: true)   : [],
        params.library ? file(params.library, checkIfExists: true) : [],
        quant
    ]))
}
""", 'process_high')

# ---------------------------------------------------------------- fastqscreen/buildfromindex
spec('fastqscreen/buildfromindex', 'FASTQSCREEN_BUILDFROMINDEX', 'fastqscreen buildfromindex', 'Build a FastQ Screen configuration from a directory of per-genome Bowtie2 index directories.',
     {'indexes': P('input', 'path', 'Parent directory holding one Bowtie2 index directory per genome; directory names become genome names.', True),
      'outdir': OUTDIR},
     {'indexes': None, 'outdir': 'results'},
     """workflow {
    dirs_ch = Channel.fromPath(params.indexes, checkIfExists: true)
        .map { d -> d.listFiles().findAll { it.isDirectory() }.toList() }

    FASTQSCREEN_BUILDFROMINDEX(dirs_ch.map { ds -> ds*.name }, dirs_ch)
}
""", 'process_single', args=False)

# ---------------------------------------------------------------- folddisco/query
spec('folddisco/query', 'FOLDDISCO_QUERY', 'folddisco query', 'Query a FoldDisco index with residue positions from a PDB structure.',
     {'pdb': P('input', 'file', 'Query PDB structure.', True),
      'query': P('option', 'string', 'Comma-separated residue queries, e.g. A13,A95,A167.'),
      'index': P('input', 'path', 'FoldDisco index directory (from folddisco index).', True),
      'structures': P('input', 'path', 'Directory of PDB structures the index was built from.', True),
      'outdir': OUTDIR},
     {'pdb': TD + 'proteomics/pdb/1tim.pdb', 'query': 'A13,A95,A167',
      'index': None, 'structures': None, 'outdir': 'results'},
     """workflow {
    structures_ch = Channel.fromPath("${params.structures}/*").collect()

    FOLDDISCO_QUERY(
        Channel.value([[id: 'query'], file(params.pdb, checkIfExists: true), params.query]),
        structures_ch.map { s -> [[id: 'index'], file(params.index, checkIfExists: true), s] }
    )
}
""", 'process_low')

# ---------------------------------------------------------------- force/tileextent
spec('force/tileextent', 'FORCE_TILEEXTENT', 'force tileextent', 'Compute the FORCE datacube tile extent for an area of interest.',
     {'aoi': P('input', 'file', 'Area-of-interest vector file (e.g. GeoPackage).', True),
      'datacube_definition': P('input', 'file', 'FORCE datacube-definition.prj file.', True),
      'shapefile_dbf': P('input', 'file', 'Optional shapefile .dbf (if the AOI is a shapefile).'),
      'shapefile_prj': P('input', 'file', 'Optional shapefile .prj.'),
      'shapefile_shx': P('input', 'file', 'Optional shapefile .shx.'),
      'outdir': OUTDIR},
     {'aoi': TD + 'earth_sciences/vector/gpkg/crete.gpkg',
      'datacube_definition': TD + 'earth_sciences/datacubes/datacube-definition.prj',
      'shapefile_dbf': None, 'shapefile_prj': None, 'shapefile_shx': None,
      'outdir': 'results'},
     """workflow {
    def shp = { p -> p ? file(p, checkIfExists: true) : [] }

    FORCE_TILEEXTENT(
        file(params.aoi, checkIfExists: true),
        file(params.datacube_definition, checkIfExists: true),
        shp(params.shapefile_dbf), shp(params.shapefile_prj), shp(params.shapefile_shx)
    )
}
""", 'process_single', args=False)

# ---------------------------------------------------------------- fqtk
spec('fqtk', 'FQTK', 'fqtk demux', 'Demultiplex paired FASTQs with fqtk using a sample sheet and per-file read structures.',
     {'sample_sheet': P('input', 'file', 'Sample metadata TSV for fqtk.', True),
      'fastq_folder': P('input', 'path', 'Directory holding the FASTQ files named in fastq_names.', True),
      'fastq_names': P('option', 'string', 'Comma-separated FASTQ file names (basenames inside fastq_folder).'),
      'read_structures': P('option', 'string', 'Comma-separated read structures matching fastq_names, e.g. 8B,8B,150T,150T.'),
      'outdir': OUTDIR},
     {'sample_sheet': 'https://github.com/nf-core/test-datasets/raw/demultiplex/testdata/sim-data/fqtk_sample_metadata_subset.tsv',
      'fastq_folder': None,
      'fastq_names': 'out_L001_I1_001.fastq.gz,out_L001_I2_001.fastq.gz,out_L001_R1_001.fastq.gz,out_L001_R2_001.fastq.gz',
      'read_structures': '8B,8B,150T,150T',
      'outdir': 'results'},
     """workflow {
    def names   = params.fastq_names.tokenize(',')*.trim()
    def structs = params.read_structures.tokenize(',')*.trim()
    if (names.size() != structs.size()) {
        error "fastq_names (${names.size()}) and read_structures (${structs.size()}) must have the same number of comma-separated entries."
    }
    def pairs = [names, structs].transpose()

    FQTK(Channel.value([
        [id: 'fqtk'],
        file(params.sample_sheet, checkIfExists: true),
        file(params.fastq_folder, checkIfExists: true),
        pairs
    ]))
}
""", 'process_high')

# ---------------------------------------------------------------- gedi/price
spec('gedi/price', 'GEDI_PRICE', 'gedi price', 'Identify translated ORFs from Ribo-seq BAMs with GEDI PRICE.',
     {'bams': P('input', 'path', 'Directory of Ribo-seq BAMs (*.bam) with their .bam.bai indexes alongside.', True),
      'index': P('input', 'path', 'GEDI genome index directory (from gedi indexgenome).', True),
      'outdir': OUTDIR},
     {'bams': None, 'index': None, 'outdir': 'results'},
     """workflow {
    bams_ch = Channel.fromPath("${params.bams}/*.bam").collect()
    bais_ch = Channel.fromPath("${params.bams}/*.bam.bai").collect()

    GEDI_PRICE(
        bams_ch.combine(bais_ch).map { b, bi -> [[id: 'price', single_end: true], b, bi] },
        Channel.value([[id: 'index'], file(params.index, checkIfExists: true)])
    )
}
""", 'process_medium')

# ---------------------------------------------------------------- interproscan
spec('interproscan', 'INTERPROSCAN', 'interproscan', 'Annotate protein sequences with InterProScan (requires the InterProScan data directory).',
     {'fasta': P('input', 'file', 'Protein FASTA to annotate.', True),
      'interproscan_database': P('input', 'path', 'InterProScan data directory (staged as data/).', True),
      'outdir': OUTDIR},
     {'fasta': TD + 'proteomics/interproscan/l_arginase.faa',
      'interproscan_database': None, 'outdir': 'results'},
     """workflow {
    INTERPROSCAN(
        Channel.value([[id: 'interproscan'], file(params.fasta, checkIfExists: true)]),
        file(params.interproscan_database, checkIfExists: true)
    )
}
""", 'process_single')

# ---------------------------------------------------------------- kaiju/mkfmi
spec('kaiju/mkfmi', 'KAIJU_MKFMI', 'kaiju mkfmi', 'Build a Kaiju FM-index from a protein FASTA and NCBI taxonomy dumps.',
     {'fasta': P('input', 'file', 'Protein FASTA to index.', True),
      'nodes_dmp': P('input', 'file', 'NCBI taxonomy nodes.dmp.', True),
      'names_dmp': P('input', 'file', 'NCBI taxonomy names.dmp.', True),
      'keep_intermediate': P('option', 'boolean', 'Keep intermediate BWT files.'),
      'outdir': OUTDIR},
     {'fasta': TD + 'genomics/sarscov2/genome/proteome.fasta',
      'nodes_dmp': TD + 'delete_me/kaiju/nodes.dmp',
      'names_dmp': TD + 'delete_me/kaiju/names.dmp',
      'keep_intermediate': False, 'outdir': 'results'},
     """workflow {
    KAIJU_MKFMI(
        Channel.value([[id: 'kaiju', single_end: true], file(params.fasta, checkIfExists: true)]),
        file(params.nodes_dmp, checkIfExists: true),
        file(params.names_dmp, checkIfExists: true),
        params.keep_intermediate
    )
}
""", 'process_high')

# ---------------------------------------------------------------- krakenuniq/build
spec('krakenuniq/build', 'KRAKENUNIQ_BUILD', 'krakenuniq build', 'Build a KrakenUniq database from custom library/taxonomy directories and a seqid2taxid map.',
     {'library': P('input', 'path', 'Directory of library FASTA files (contents staged under library/).', True),
      'taxonomy': P('input', 'path', 'Directory of taxonomy dump files (contents staged under taxonomy/).', True),
      'seqid2taxid': P('input', 'file', 'seqid2taxid.map file.', True),
      'keep_intermediate': P('option', 'boolean', 'Keep intermediate database files.'),
      'outdir': OUTDIR},
     {'library': None, 'taxonomy': None,
      'seqid2taxid': TD + 'genomics/sarscov2/metagenome/seqid2taxid.map',
      'keep_intermediate': False, 'outdir': 'results'},
     """workflow {
    lib_ch = Channel.fromPath("${params.library}/*").collect()
    tax_ch = Channel.fromPath("${params.taxonomy}/*").collect()

    KRAKENUNIQ_BUILD(
        lib_ch.combine(tax_ch).map { lib, tax ->
            [[id: 'krakenuniq'], lib, tax, file(params.seqid2taxid, checkIfExists: true)]
        },
        params.keep_intermediate
    )
}
""", 'process_medium')

# ---------------------------------------------------------------- krona/ktimporttaxonomy
spec('krona/ktimporttaxonomy', 'KRONA_KTIMPORTTAXONOMY', 'krona ktimporttaxonomy', 'Convert a Kraken report to a Krona HTML plot using a Krona taxonomy table.',
     {'report': P('input', 'file', 'Kraken/Kraken2 report file.', True),
      'taxonomy': P('input', 'file', 'Krona taxonomy.tab table.', True),
      'outdir': OUTDIR},
     {'report': TD + 'genomics/sarscov2/metagenome/test_1.kraken2.report.txt',
      'taxonomy': TD + 'genomics/sarscov2/metagenome/krona_taxonomy.tab',
      'outdir': 'results'},
     """workflow {
    KRONA_KTIMPORTTAXONOMY(
        Channel.value([[id: 'krona'], file(params.report, checkIfExists: true)]),
        file(params.taxonomy, checkIfExists: true)
    )
}
""", 'process_single')

# ---------------------------------------------------------------- merqury/hapmers
spec('merqury/hapmers', 'MERQURY_HAPMERS', 'merqury hapmers', 'Find parent-specific hap-mers from child/maternal/paternal meryl databases.',
     {'child_meryl': P('input', 'path', 'Child meryl database directory.', True),
      'maternal_meryl': P('input', 'path', 'Maternal meryl database directory.', True),
      'paternal_meryl': P('input', 'path', 'Paternal meryl database directory.', True),
      'outdir': OUTDIR},
     {'child_meryl': None, 'maternal_meryl': None, 'paternal_meryl': None, 'outdir': 'results'},
     """workflow {
    MERQURY_HAPMERS(
        Channel.value([[id: 'hapmers'], file(params.child_meryl, checkIfExists: true)]),
        file(params.maternal_meryl, checkIfExists: true),
        file(params.paternal_meryl, checkIfExists: true)
    )
}
""", 'process_single')

# ---------------------------------------------------------------- metabuli/build
spec('metabuli/build', 'METABULI_BUILD', 'metabuli build', 'Build a Metabuli taxonomy database from sequences and NCBI taxonomy dumps.',
     {'fasta': P('input', 'file', 'Sequence FASTA to build the database from.', True),
      'taxonomy_names': P('input', 'file', 'NCBI taxonomy names.dmp.', True),
      'taxonomy_nodes': P('input', 'file', 'NCBI taxonomy nodes.dmp.', True),
      'taxonomy_merged': P('input', 'file', 'Optional NCBI taxonomy merged.dmp.'),
      'accession2taxid': P('input', 'file', 'accession2taxid mapping file (e.g. nucl_gb.accession2taxid).', True),
      'cds_info': P('input', 'file', 'Optional CDS info file.'),
      'outdir': OUTDIR},
     {'fasta': TD + 'genomics/sarscov2/genome/genome.fasta',
      'taxonomy_names': TD + 'genomics/prokaryotes/metagenome/taxonomy/taxdmp/names.dmp',
      'taxonomy_nodes': TD + 'genomics/prokaryotes/metagenome/taxonomy/taxdmp/nodes.dmp',
      'taxonomy_merged': None,
      'accession2taxid': TD + 'genomics/prokaryotes/metagenome/taxonomy/accession2taxid/nucl_gb.accession2taxid',
      'cds_info': None, 'outdir': 'results'},
     """workflow {
    def opt = { p -> p ? file(p, checkIfExists: true) : [] }

    METABULI_BUILD(
        Channel.value([[id: 'metabuli'], [file(params.fasta, checkIfExists: true)]]),
        file(params.taxonomy_names, checkIfExists: true),
        file(params.taxonomy_nodes, checkIfExists: true),
        opt(params.taxonomy_merged),
        [file(params.accession2taxid, checkIfExists: true)],
        opt(params.cds_info)
    )
}
""", 'process_medium')

# ---------------------------------------------------------------- mmseqs/createtsv
spec('mmseqs/createtsv', 'MMSEQS_CREATETSV', 'mmseqs createtsv', 'Convert an MMseqs2 result database to a TSV alignment table.',
     {'db_result': P('input', 'path', 'MMseqs2 result database directory (from cluster/search/taxonomy).', True),
      'db_query': P('input', 'path', 'Optional query database directory (from mmseqs createdb).'),
      'db_target': P('input', 'path', 'Target database directory (from mmseqs createdb).', True),
      'outdir': OUTDIR},
     {'db_result': None, 'db_query': None, 'db_target': None, 'outdir': 'results'},
     """workflow {
    MMSEQS_CREATETSV(
        Channel.value([[id: 'createtsv'], file(params.db_result, checkIfExists: true)]),
        Channel.value([[id: 'query'], params.db_query ? file(params.db_query, checkIfExists: true) : []]),
        Channel.value([[id: 'target'], file(params.db_target, checkIfExists: true)])
    )
}
""", 'process_single')

# ---------------------------------------------------------------- multiqcsav
spec('multiqcsav', 'MULTIQCSAV', 'multiqcsav', 'Generate a MultiQC report from an Illumina sequencing run (SAV InterOp data).',
     {'run_dir': P('input', 'path', 'Illumina run directory containing RunInfo.xml and InterOp/*.bin.', True),
      'extra_multiqc_files': P('input', 'path', 'Optional directory of extra files for MultiQC (e.g. FastQC zips).'),
      'multiqc_config': P('input', 'file', 'Optional MultiQC config YAML.'),
      'multiqc_logo': P('input', 'file', 'Optional custom report logo.'),
      'replace_names': P('input', 'file', 'Optional MultiQC sample-name replacement TSV.'),
      'sample_names': P('input', 'file', 'Optional MultiQC sample-name cleaning file.'),
      'outdir': OUTDIR},
     {'run_dir': None, 'extra_multiqc_files': None,
      'multiqc_config': 'https://raw.githubusercontent.com/nf-core/seqinspector/1.0.0/assets/multiqc_config.yml',
      'multiqc_logo': None, 'replace_names': None, 'sample_names': None, 'outdir': 'results'},
     """workflow {
    def run = file(params.run_dir, checkIfExists: true)
    def xml = []
    def interop = []
    run.eachFileRecurse { f ->
        if (f.fileName.toString() == 'RunInfo.xml') xml << f
        if (f.parent?.fileName?.toString() == 'InterOp' && f.fileName.toString().endsWith('.bin')) interop << f
    }
    def extra = params.extra_multiqc_files
        ? file(params.extra_multiqc_files, checkIfExists: true).listFiles().findAll { it.isFile() }.toList()
        : []
    def opt = { p -> p ? file(p, checkIfExists: true) : [] }

    MULTIQCSAV(Channel.value([
        [id: 'sav'], xml, interop, extra,
        opt(params.multiqc_config), opt(params.multiqc_logo),
        opt(params.replace_names), opt(params.sample_names)
    ]))
}
""", 'process_single')

# ---------------------------------------------------------------- ribotish/predict
spec('ribotish/predict', 'RIBOTISH_PREDICT', 'ribotish predict', 'Predict translated ORFs from Ribo-seq (and optional TI-seq) BAMs with RiboTISH.',
     {'bam_ribo': P('input', 'file', 'Ribo-seq aligned BAM.', True),
      'bai_ribo': P('input', 'file', 'Ribo-seq BAM index.', True),
      'bam_ti': P('input', 'file', 'Optional TI-seq aligned BAM.'),
      'bai_ti': P('input', 'file', 'Optional TI-seq BAM index.'),
      'fasta': P('input', 'file', 'Genome FASTA (uncompressed).', True),
      'gtf': P('input', 'file', 'Gene annotation GTF.', True),
      'candidate_orfs': P('input', 'file', 'Optional candidate ORFs file from a previous ribotish run.'),
      'para_ribo': P('input', 'file', 'Optional Ribo-seq parameter file.'),
      'para_ti': P('input', 'file', 'Optional TI-seq parameter file.'),
      'reference_gtf': P('input', 'file', 'Optional secondary reference GTF.'),
      'outdir': OUTDIR},
     {'bam_ribo': TD + 'genomics/homo_sapiens/riboseq_expression/aligned_reads/SRX11780888_chr20.bam',
      'bai_ribo': TD + 'genomics/homo_sapiens/riboseq_expression/aligned_reads/SRX11780888_chr20.bam.bai',
      'bam_ti': None, 'bai_ti': None, 'fasta': None,
      'gtf': TD + 'genomics/homo_sapiens/riboseq_expression/Homo_sapiens.GRCh38.111_chr20.gtf',
      'candidate_orfs': None, 'para_ribo': None, 'para_ti': None, 'reference_gtf': None,
      'outdir': 'results'},
     """workflow {
    def opt2 = { a, b -> (a && b) ? [file(a, checkIfExists: true), file(b, checkIfExists: true)] : [[], []] }
    def opt1 = { p -> p ? [file(p, checkIfExists: true)] : [[]] }

    RIBOTISH_PREDICT(
        Channel.value([[id: 'ribo', single_end: true, strandedness: 'forward'],
                       file(params.bam_ribo, checkIfExists: true), file(params.bai_ribo, checkIfExists: true)]),
        Channel.value([[id: 'ti']] + opt2(params.bam_ti, params.bai_ti)),
        Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.gtf, checkIfExists: true)]),
        Channel.value([[id: 'orfs']] + opt1(params.candidate_orfs)),
        Channel.value([[id: 'para_ribo']] + opt1(params.para_ribo)),
        Channel.value([[id: 'para_ti']] + opt1(params.para_ti)),
        Channel.value([[id: 'secondary']] + opt1(params.reference_gtf))
    )
}
""", 'process_medium')

# ---------------------------------------------------------------- riker/multi
spec('riker/multi', 'RIKER_MULTI', 'riker multi', 'Run RIKER multi-tool QC (alignment, basic, isize, ...) on a BAM.',
     {'bam': P('input', 'file', 'Coordinate-sorted BAM.', True),
      'bai': P('input', 'file', 'BAM index.', True),
      'error_vcf': P('input', 'file', 'Optional error-model VCF.'),
      'error_vcf_idx': P('input', 'file', 'Optional error-model VCF index.'),
      'error_intervals': P('input', 'file', 'Optional error-model intervals.'),
      'gcbias_exclude_intervals': P('input', 'file', 'Optional GC-bias exclusion intervals.'),
      'hybcap_baits': P('input', 'path', 'Optional hybrid-capture baits directory.'),
      'hybcap_targets': P('input', 'path', 'Optional hybrid-capture targets directory.'),
      'rna_gene_model': P('input', 'file', 'Optional RNA gene model (refFlat).'),
      'rna_ribosomal_intervals': P('input', 'file', 'Optional ribosomal interval list.'),
      'wgs_intervals': P('input', 'file', 'Optional WGS interval list.'),
      'fasta': P('input', 'file', 'Optional reference FASTA.'),
      'fai': P('input', 'file', 'Optional reference FASTA index.'),
      'args': P('option', 'string', 'Extra riker arguments, e.g. --tools alignment basic isize.'),
      'outdir': OUTDIR},
     {'bam': TD + 'genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam',
      'bai': TD + 'genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam.bai',
      'error_vcf': None, 'error_vcf_idx': None, 'error_intervals': None,
      'gcbias_exclude_intervals': None, 'hybcap_baits': None, 'hybcap_targets': None,
      'rna_gene_model': None, 'rna_ribosomal_intervals': None, 'wgs_intervals': None,
      'fasta': None, 'fai': None, 'args': '--tools alignment basic isize', 'outdir': 'results'},
     """workflow {
    def opt = { p -> p ? file(p, checkIfExists: true) : [] }
    def dirlist = { p -> p ? file(p, checkIfExists: true).listFiles().findAll { it.isFile() }.toList() : [] }

    RIKER_MULTI(
        Channel.value([
            [id: 'riker', single_end: false],
            file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true),
            opt(params.error_vcf), opt(params.error_vcf_idx), opt(params.error_intervals),
            opt(params.gcbias_exclude_intervals),
            dirlist(params.hybcap_baits), dirlist(params.hybcap_targets),
            opt(params.rna_gene_model), opt(params.rna_ribosomal_intervals), opt(params.wgs_intervals)
        ]),
        Channel.value([[id: 'reference'], opt(params.fasta), opt(params.fai)])
    )
}
""", 'process_medium')

# ---------------------------------------------------------------- samplesheetparser/diff
spec('samplesheetparser/diff', 'SAMPLESHEETPARSER_DIFF', 'samplesheetparser diff', 'Diff two Illumina sample sheets and report the changes.',
     {'old_sheet': P('input', 'file', 'Original sample sheet CSV.', True),
      'new_sheet': P('input', 'file', 'Updated sample sheet CSV.', True),
      'outdir': OUTDIR},
     {'old_sheet': TD + 'genomics/homo_sapiens/illumina/bcl/flowcell_samplesheet.v1.csv',
      'new_sheet': TD + 'genomics/homo_sapiens/illumina/bcl/flowcell_samplesheet.v1.csv',
      'outdir': 'results'},
     """workflow {
    SAMPLESHEETPARSER_DIFF(Channel.value([
        [id: 'diff'],
        file(params.old_sheet, checkIfExists: true),
        file(params.new_sheet, checkIfExists: true)
    ]))
}
""", 'process_single')

# ---------------------------------------------------------------- samplesheetparser/merge
spec('samplesheetparser/merge', 'SAMPLESHEETPARSER_MERGE', 'samplesheetparser merge', 'Merge multiple Illumina sample sheets into one.',
     {'samplesheets': P('input', 'path', 'Directory of sample sheet CSVs (*.csv).', True),
      'target_version': P('option', 'string', 'Sample sheet format version to emit.', enum=['v1', 'v2']),
      'outdir': OUTDIR},
     {'samplesheets': None, 'target_version': 'v1', 'outdir': 'results'},
     """workflow {
    sheets_ch = Channel.fromPath("${params.samplesheets}/*.csv").collect()

    SAMPLESHEETPARSER_MERGE(sheets_ch.map { s -> [[id: 'merge'], s] }, params.target_version)
}
""", 'process_single')

# ---------------------------------------------------------------- sourcepredict
spec('sourcepredict', 'SOURCEPREDICT', 'sourcepredict', 'Predict contamination sources of a metagenomic sample with SourcePredict.',
     {'kraken_parse': P('input', 'file', 'Parsed Kraken report of the sink sample (taxid CSV).', True),
      'sources': P('input', 'file', 'Source samples taxid CSV.', True),
      'labels': P('input', 'file', 'Source labels CSV.', True),
      'taxa_sqlite': P('input', 'file', 'NCBI taxa.sqlite database (decompressed).', True),
      'taxa_sqlite_traverse_pkl': P('input', 'file', 'Precomputed taxonomy traversal pickle.', True),
      'outdir': OUTDIR},
     {'kraken_parse': TD + 'genomics/prokaryotes/metagenome/taxonomy/misc/sink_taxid.csv',
      'sources': TD + 'genomics/prokaryotes/metagenome/taxonomy/misc/sources_taxid.csv',
      'labels': TD + 'genomics/prokaryotes/metagenome/taxonomy/misc/sources_labels.csv',
      'taxa_sqlite': None,
      'taxa_sqlite_traverse_pkl': TD + 'genomics/prokaryotes/metagenome/taxonomy/misc/taxa_sqlite_traverse.pkl',
      'outdir': 'results'},
     """workflow {
    SOURCEPREDICT(
        Channel.value([[id: 'sourcepredict'], file(params.kraken_parse, checkIfExists: true)]),
        file(params.sources, checkIfExists: true),
        file(params.labels, checkIfExists: true),
        file(params.taxa_sqlite, checkIfExists: true),
        file(params.taxa_sqlite_traverse_pkl, checkIfExists: true)
    )
}
""", 'process_medium')

# ---------------------------------------------------------------- stitch
spec('stitch', 'STITCH', 'STITCH', 'Impute low-coverage sequencing data with STITCH (BAM or CRAM inputs).',
     {'crams': P('input', 'path', 'Directory of per-sample CRAM (or BAM) files.', True),
      'crais': P('input', 'path', 'Directory of the matching .crai (or .bai) indexes.', True),
      'posfile': P('input', 'file', 'Position file (known variant sites).', True),
      'input': P('input', 'file', 'Optional pre-generated STITCH input file.'),
      'genetic_map': P('input', 'file', 'Optional genetic map file.'),
      'rdata': P('input', 'path', 'Optional RData directory from a previous STITCH run.'),
      'samplename': P('input', 'file', 'Optional sample names file.'),
      'chromosome_name': P('option', 'string', 'Chromosome to impute.'),
      'start': P('option', 'integer', 'Region start (optional).'),
      'end': P('option', 'integer', 'Region end (optional).'),
      'buffer': P('option', 'integer', 'Region buffer.'),
      'k': P('option', 'integer', 'Number of ancestral haplotypes (K).'),
      'n_gen': P('option', 'integer', 'Number of generations since founding (nGen).'),
      'fasta': P('input', 'file', 'Reference FASTA (required for CRAM).', True),
      'fai': P('input', 'file', 'Reference FASTA index.', True),
      'seed': P('option', 'integer', 'Random seed for reproducibility.'),
      'outdir': OUTDIR},
     {'crams': None, 'crais': None,
      'posfile': TD + 'genomics/homo_sapiens/genome/chr21/sequence/dbsnp_138.hg38.first_10_biallelic_sites.tsv',
      'input': TD + 'genomics/homo_sapiens/genome/chr21/sequence/genome.GRCh38.stitch.chr21.txt.gz',
      'genetic_map': None, 'rdata': None, 'samplename': None,
      'chromosome_name': 'chr21', 'start': 8522360, 'end': 8522670, 'buffer': 0,
      'k': 2, 'n_gen': 1,
      'fasta': TD + 'genomics/homo_sapiens/genome/chr21/sequence/genome.fasta',
      'fai': TD + 'genomics/homo_sapiens/genome/chr21/sequence/genome.fasta.fai',
      'seed': 1, 'outdir': 'results'},
     """workflow {
    ch_reads = Channel.fromPath("${params.crams}/*").collect()
        .map { fs -> fs.findAll { it.name.endsWith('.cram') || it.name.endsWith('.bam') } }
    ch_idx = Channel.fromPath("${params.crais}/*").collect()
    ch_list = ch_reads.map { fs -> fs*.name.join('\\n') }.collectFile(name: 'cramlist.txt', newLine: true)
    def opt = { p -> p ? file(p, checkIfExists: true) : [] }

    ch_main = ch_reads.combine(ch_idx).combine(ch_list).map { reads, idx, list ->
        [[id: 'stitch'], reads, idx, list,
         opt(params.samplename),
         file(params.posfile, checkIfExists: true),
         opt(params.input),
         opt(params.genetic_map),
         opt(params.rdata),
         params.chromosome_name, params.start, params.end, params.buffer, params.k, params.n_gen]
    }

    STITCH(ch_main,
           Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)]),
           params.seed)
}
""", 'process_medium')

# ---------------------------------------------------------------- td2/predict
spec('td2/predict', 'TD2_PREDICT', 'td2 predict', 'Predict translated ORFs with TD2 using a genome FASTA and a LongOrfs ORF directory.',
     {'fasta': P('input', 'file', 'Genome FASTA.', True),
      'orfs_dir': P('input', 'path', 'ORF directory from td2 longorfs.', True),
      'outdir': OUTDIR},
     {'fasta': TD + 'genomics/homo_sapiens/genome/genome.fasta',
      'orfs_dir': None, 'outdir': 'results'},
     """workflow {
    TD2_PREDICT(Channel.value([
        [id: 'td2'],
        file(params.fasta, checkIfExists: true),
        file(params.orfs_dir, checkIfExists: true)
    ]))
}
""", 'process_high')

# ---------------------------------------------------------------- vamb/bin
spec('vamb/bin', 'VAMB_BIN', 'vamb bin', 'Bin metagenomic contigs with VAMB using co-abundance BAMs.',
     {'assembly': P('input', 'file', 'Assembly contigs FASTA (optionally gzipped).', True),
      'abundance_tsv': P('input', 'file', 'Optional precomputed abundance TSV (instead of BAMs).'),
      'bams': P('input', 'path', 'Directory of per-sample sorted BAMs (*.bam).'),
      'taxonomy': P('input', 'file', 'Optional taxonomy file for taxonomic binning.'),
      'outdir': OUTDIR},
     {'assembly': TD + 'genomics/prokaryotes/metagenome/binning/multisample/contigs.fasta.gz',
      'abundance_tsv': None, 'bams': None, 'taxonomy': None, 'outdir': 'results'},
     """workflow {
    def bams = params.bams
        ? file(params.bams, checkIfExists: true).listFiles().findAll { it.isFile() && it.name.endsWith('.bam') }.toList()
        : []
    def opt = { p -> p ? file(p, checkIfExists: true) : [] }

    VAMB_BIN(Channel.value([
        [id: 'vamb'],
        file(params.assembly, checkIfExists: true),
        opt(params.abundance_tsv),
        bams,
        opt(params.taxonomy)
    ]))
}
""", 'process_high')


def main():
    assert len(SPECS) == 27, len(SPECS)
    for s in SPECS:
        rel, proc = s['rel'], s['proc']
        if s['args'] and 'args' not in s['params']:
            tool = s['name'].split()[0].lower()
            items = [(k, v) for k, v in s['params'].items() if k != 'outdir']
            items.append(('args', P('option', 'string', f'Extra command-line arguments passed to {tool}.')))
            items.append(('outdir', s['params']['outdir']))
            s['params'] = dict(items)
            s['defaults'] = {**{k: v for k, v in s['defaults'].items() if k != 'outdir'},
                             'args': '', 'outdir': s['defaults']['outdir']}
        wdir = os.path.join(DST, rel, 'wrapper')
        if os.path.exists(os.path.join(wdir, 'wrapper.yaml')):
            print('skip (exists)', rel)
            continue
        os.makedirs(wdir, exist_ok=True)
        # params declarations in main.nf
        decl = ''.join(f'params.{k} = null\n' for k in s['params'])
        note = (s['note'] + '\n') if s['note'] else ''
        main_nf = HEADER % {'proc': proc} + note + decl + '\n' + s['body']
        open(os.path.join(wdir, 'main.nf'), 'w').write(main_nf)
        # wrapper.yaml
        y = {'id': 'nf-core/modules/' + rel.replace('/', '-').replace('_', '-'),
             'name': s['name'], 'summary': s['summary'],
             'params': s['params'],
             'outputs': {'output': {'type': 'directory',
                                    'path': '${outdir}/' + proc.lower(),
                                    'primary': True}}}
        with open(os.path.join(wdir, 'wrapper.yaml'), 'w') as f:
            f.write('# Agent-facing wrapper adapter for the vendored ' + rel + ' module at\n')
            f.write('# ../main.nf. Follows the minimal contract in\n')
            f.write('# docs/design/phi-wrapper-agent-composition-design.md section 3.\n')
            f.write('# Manually wired: exotic staging pattern.\n')
            yaml.safe_dump(y, f, sort_keys=False, allow_unicode=True, width=100)
        # params.json
        defaults = dict(s['defaults'])
        assert set(defaults) == set(s['params']), (rel, set(s['params']) - set(defaults), set(defaults) - set(s['params']))
        open(os.path.join(wdir, 'params.json'), 'w').write(json.dumps(defaults, indent=2) + '\n')
        # nextflow.config
        cpus, mem, time = RESOURCES[s['label']]
        extargs = "        ext.args = { params.args ?: '' }\n" if s['args'] else ''
        cfg = CONFIG % {'proc': proc, 'cpus': cpus, 'mem': mem, 'time': time,
                        'outdir_name': proc.lower(), 'extargs': extargs}
        open(os.path.join(wdir, 'nextflow.config'), 'w').write(cfg)
        print('wired', rel)

if __name__ == '__main__':
    main()