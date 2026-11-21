#!/usr/bin/env python3
"""Port the Snakemake-only tools (no nf-core equivalent) to native Nextflow local modules + phi wrappers.

Source: /Users/dengxsh/Documents/WorkSpace/VSCode/Workflow/wrappers (MIT, dengxsh).
metawrap is skipped: its Snakemake wrapper dirs are empty placeholders.
"""
import json, os, re
import yaml

DST = 'resources/wrappers/modules/local'
QUAY = 'quay.io/biocontainers/'

RESOURCES = {
    'process_single': (1, '2.GB', '1.h'),
    'process_low': (2, '4.GB', '2.h'),
    'process_medium': (4, '8.GB', '4.h'),
    'process_high': (8, '16.GB', '8.h'),
    'process_long': (2, '4.GB', '20.h'),
    'process_high_memory': (8, '64.GB', '8.h'),
}

def container(tag):
    return ("    container \"${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?\n"
            f"        'docker://{QUAY}{tag}' :\n"
            f"        '{QUAY}{tag}' }}\"\n")

def versions(tool, ver):
    return f'''
    cat <<-END_VERSIONS > versions.yml
    "${{task.process}}":
        {tool}: {ver}
    END_VERSIONS
    '''

MODULES = {}

def mod(rel, proc, label, tag, inputs, outputs, script, stub, env_deps, meta, wrapper):
    MODULES[rel] = dict(proc=proc, label=label, tag=tag, inputs=inputs, outputs=outputs,
                        script=script, stub=stub, env_deps=env_deps, meta=meta, wrapper=wrapper)

def P(kind, typ, desc, required=False, enum=None):
    d = {'kind': kind, 'type': typ, 'required': required, 'description': desc}
    if enum: d['enum'] = enum
    return d

OUTDIR = P('output', 'path', 'Output directory for the results.', True)
ARGS = P('option', 'string', 'Extra command-line arguments passed to the tool.')

def meta_yml(name, desc, keywords, tool, tool_desc, homepage, inputs, outputs):
    return {
        'name': name,
        'description': desc,
        'keywords': keywords,
        'tools': [{tool: {'description': tool_desc, 'homepage': homepage}}],
        'input': inputs,
        'output': outputs,
        'authors': ['@dengxsh'],
        'maintainers': ['@dengxsh'],
    }

# ---------------------------------------------------------------- 3d-dna
mod('3d-dna', 'THREEDDNA', 'process_high', '3d-dna:201008--hdfd78af_0',
    """    tuple val(meta), path(fasta), path(mnd)
""",
    """    tuple val(meta), path("*.Final.fasta")           , emit: final_fasta, optional: true
    tuple val(meta), path("*.final.assembly.fasta")  , emit: assembly, optional: true
    tuple val(meta), path("*.hic")                   , emit: hic, optional: true
    path "versions.yml"                              , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    """ + '''
    """
    3d-dna \\\\
        ${args} \\\\
        ${fasta} \\\\
        ${mnd}
    """ + versions('3d-dna', '201008') + ''',
    """    touch ${meta.id}.Final.fasta
    touch ${meta.id}.final.assembly.fasta
    touch ${meta.id}.hic
    """ + versions('3d-dna', '201008'),
    ['3d-dna=201008'],
    meta_yml('3d-dna', 'Scaffold a draft genome assembly with 3D-DNA using Hi-C contact maps.',
             ['scaffolding', 'hi-c', 'genome assembly'],
             '3d-dna', 'The 3D-DNA scaffolding pipeline.', 'https://github.com/aidenlab/3d-dna',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'fasta': {'type': 'file', 'description': 'Draft assembly FASTA', 'pattern': '*.{fasta,fa,fna}'}},
              {'mnd': {'type': 'file', 'description': 'Merged nodup Hi-C contacts (merged_nodups.txt from Juicer)', 'pattern': '*.{txt,mnd}'}}],
             [{'final_fasta': {'type': 'file', 'description': 'Final scaffolded FASTA'}},
              {'assembly': {'type': 'file', 'description': 'Final assembly FASTA'}}]),
    dict(name='3D-DNA', summary='Scaffold a draft genome with 3D-DNA using Hi-C merged_nodups contacts.',
         params={'fasta': P('input', 'file', 'Draft assembly FASTA.', True),
                 'mnd': P('input', 'file', 'merged_nodups.txt Hi-C contact file (from Juicer).', True),
                 'mode': P('option', 'string', 'Assembly mode.', enum=['haploid', 'diploid']),
                 'input_size': P('option', 'integer', 'Input contig size threshold (-i).'),
                 'rounds': P('option', 'integer', 'Number of misjoin-correction rounds (-r).'),
                 'stage': P('option', 'string', 'Stage to start from (-s), e.g. seal, finalize.'),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'fasta': None, 'mnd': None, 'mode': 'haploid', 'input_size': 15000,
                   'rounds': 2, 'stage': None, 'args': '', 'outdir': 'results'},
         body="""workflow {
    THREEDDNA(Channel.value([
        [id: 'threeddna'],
        file(params.fasta, checkIfExists: true),
        file(params.mnd, checkIfExists: true)
    ]))
}
""",
         extargs="""        ext.args = { "--mode ${params.mode} -i ${params.input_size} -r ${params.rounds}" +
            (params.stage ? " -s ${params.stage}" : '') + " ${params.args ?: ''}" }
"""))

# ---------------------------------------------------------------- dramv
mod('dramv', 'DRAMV', 'process_high', 'dram:1.5.0--pyhdfd78af_0',
    """    tuple val(meta), path(fasta), path(affi)
    path(config)
    val(distill)
""",
    """    tuple val(meta), path("annotate/annotations.tsv"), emit: annotations
    tuple val(meta), path("distill")                 , emit: distill, optional: true
    path "versions.yml"                              , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def config_opt = config ? "--config_loc ${config}" : ''
    """ + '''
    """
    DRAM-v.py annotate \\\\
        --threads ${task.cpus} \\\\
        -i ${fasta} \\\\
        -v ${affi} \\\\
        -o annotate \\\\
        ${config_opt} \\\\
        ${args}

    if [ "${distill}" = "true" ]; then
        DRAM-v.py distill \\\\
            -i annotate/annotations.tsv \\\\
            -o distill \\\\
            ${config_opt}
    fi
    """ + versions('dram', '1.5.0') + ''',
    """    mkdir -p annotate distill
    touch annotate/annotations.tsv
    """ + versions('dram', '1.5.0'),
    ['dram=1.5.0', 'setuptools>=59,<81', 'python>=3.9,<3.11'],
    meta_yml('dramv', 'Annotate viral contigs with DRAM-v (annotate + optional distill), from VirSorter2 for-dramv outputs.',
             ['virus', 'annotation', 'metagenomics', 'amg'],
             'dram', 'DRAM-v: distillation of viral auxiliary metabolic genes.', 'https://github.com/WrightonLabCSU/DRAM',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'fasta': {'type': 'file', 'description': 'final-viral-combined-for-dramv.fa (VirSorter2 for-dramv output)'}},
              {'affi': {'type': 'file', 'description': 'viral-affi-contigs-for-dramv.tab (keeps VirSorter2 > headers)'}},
              {'config': {'type': 'file', 'description': 'Optional DRAM export_config file with a configured amg_database'}},
              {'distill': {'type': 'boolean', 'description': 'Also run DRAM-v.py distill on the annotations'}}],
             [{'annotations': {'type': 'file', 'description': 'DRAM-v annotations.tsv'}},
              {'distill': {'type': 'directory', 'description': 'DRAM-v distill output directory'}}]),
    dict(name='DRAM-v', summary='Annotate viral contigs and AMGs with DRAM-v (annotate + distill) from VirSorter2 for-dramv inputs.',
         params={'fasta': P('input', 'file', 'for-dramv viral FASTA (final-viral-combined-for-dramv.fa).', True),
                 'affi': P('input', 'file', 'for-dramv affinity table (viral-affi-contigs-for-dramv.tab).', True),
                 'config': P('input', 'file', 'DRAM export_config file with a configured AMG database.'),
                 'distill': P('option', 'boolean', 'Also run DRAM-v.py distill.'),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'fasta': None, 'affi': None, 'config': None, 'distill': True, 'args': '', 'outdir': 'results'},
         body="""// Requires a DRAM installation with databases configured (DRAM-setup.py
// prepare_databases) and a config file whose amg_database is a real path.
workflow {
    DRAMV(
        Channel.value([[id: 'dramv'],
                       file(params.fasta, checkIfExists: true),
                       file(params.affi, checkIfExists: true)]),
        params.config ? file(params.config, checkIfExists: true) : [],
        params.distill
    )
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- graphlan/annotate
mod('graphlan/annotate', 'GRAPHLAN_ANNOTATE', 'process_single', 'graphlan:1.1.3--2',
    """    tuple val(meta), path(tree), path(annotation)
""",
    """    tuple val(meta), path("${prefix}.xml"), emit: xml
    path "versions.yml"                   , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ + '''
    """
    graphlan_annotate.py \\\\
        --annot ${annotation} \\\\
        ${args} \\\\
        ${tree} \\\\
        ${prefix}.xml
    """ + versions('graphlan', '1.1.3') + ''',
    """    touch ${meta.id}.xml
    """ + versions('graphlan', '1.1.3'),
    ['graphlan=1.1.3', 'export2graphlan=0.22'],
    meta_yml('graphlan_annotate', 'Apply a GraPhlAn annotation file to a phylogenetic tree.',
             ['visualization', 'phylogeny', 'annotation'],
             'graphlan', 'GraPhlAn: compact phylogenetic tree visualization.', 'https://github.com/biobakery/graphlan',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'tree': {'type': 'file', 'description': 'Input tree (Newick or GraPhlAn XML)'}},
              {'annotation': {'type': 'file', 'description': 'GraPhlAn annotation file'}}],
             [{'xml': {'type': 'file', 'description': 'Annotated GraPhlAn XML tree'}}]),
    dict(name='graphlan annotate', summary='Apply a GraPhlAn annotation file to a phylogenetic tree.',
         params={'tree': P('input', 'file', 'Input tree (Newick or GraPhlAn XML).', True),
                 'annotation': P('input', 'file', 'GraPhlAn annotation file.', True),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'tree': None, 'annotation': None, 'args': '', 'outdir': 'results'},
         body="""workflow {
    GRAPHLAN_ANNOTATE(Channel.value([
        [id: 'graphlan'],
        file(params.tree, checkIfExists: true),
        file(params.annotation, checkIfExists: true)
    ]))
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- graphlan/export
mod('graphlan/export', 'GRAPHLAN_EXPORT', 'process_single', 'export2graphlan:0.22--py_0',
    """    tuple val(meta), path(taxonomy)
""",
    """    tuple val(meta), path("${prefix}.tree.xml") , emit: tree
    tuple val(meta), path("${prefix}.annot.xml"), emit: annotation
    path "versions.yml"                         , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ + '''
    """
    export2graphlan.py \\\\
        -i ${taxonomy} \\\\
        --tree ${prefix}.tree.xml \\\\
        --annotation ${prefix}.annot.xml \\\\
        ${args}
    """ + versions('export2graphlan', '0.22') + ''',
    """    touch ${meta.id}.tree.xml
    touch ${meta.id}.annot.xml
    """ + versions('export2graphlan', '0.22'),
    ['graphlan=1.1.3', 'export2graphlan=0.22'],
    meta_yml('graphlan_export', 'Convert a taxonomy abundance table to a GraPhlAn tree and annotation with export2graphlan.',
             ['visualization', 'taxonomy', 'metagenomics'],
             'export2graphlan', 'Export taxonomic profiles to GraPhlAn.', 'https://github.com/biobakery/export2graphlan',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'taxonomy': {'type': 'file', 'description': 'Taxonomy abundance table (e.g. MetaPhlAn profile)'}}],
             [{'tree': {'type': 'file', 'description': 'GraPhlAn tree XML'}},
              {'annotation': {'type': 'file', 'description': 'GraPhlAn annotation XML'}}]),
    dict(name='export2graphlan', summary='Convert a taxonomy abundance table into a GraPhlAn tree + annotation.',
         params={'taxonomy': P('input', 'file', 'Taxonomy abundance table (e.g. MetaPhlAn profile).', True),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'taxonomy': None, 'args': '', 'outdir': 'results'},
         body="""workflow {
    GRAPHLAN_EXPORT(Channel.value([[id: 'graphlan'], file(params.taxonomy, checkIfExists: true)]))
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- graphlan/graphlan
mod('graphlan/graphlan', 'GRAPHLAN_GRAPHLAN', 'process_single', 'graphlan:1.1.3--2',
    """    tuple val(meta), path(xml)
    val(format)
""",
    """    tuple val(meta), path("${prefix}.${format}"), emit: image
    path "versions.yml"                         , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ + '''
    """
    graphlan.py \\\\
        ${args} \\\\
        --format ${format} \\\\
        ${xml} \\\\
        ${prefix}.${format}
    """ + versions('graphlan', '1.1.3') + ''',
    """    touch ${meta.id}.${format}
    """ + versions('graphlan', '1.1.3'),
    ['graphlan=1.1.3', 'export2graphlan=0.22'],
    meta_yml('graphlan', 'Render an annotated GraPhlAn XML tree to a circular phylogenetic figure.',
             ['visualization', 'phylogeny'],
             'graphlan', 'GraPhlAn: compact phylogenetic tree visualization.', 'https://github.com/biobakery/graphlan',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'xml': {'type': 'file', 'description': 'Annotated GraPhlAn XML tree'}},
              {'format': {'type': 'string', 'description': 'Output image format (png, svg, pdf)'}}],
             [{'image': {'type': 'file', 'description': 'Rendered tree figure'}}]),
    dict(name='graphlan', summary='Render an annotated GraPhlAn XML tree to a circular phylogenetic figure.',
         params={'xml': P('input', 'file', 'Annotated GraPhlAn XML tree.', True),
                 'format': P('option', 'string', 'Output image format.', enum=['png', 'svg', 'pdf']),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'xml': None, 'format': 'png', 'args': '', 'outdir': 'results'},
         body="""workflow {
    GRAPHLAN_GRAPHLAN(Channel.value([[id: 'graphlan'], file(params.xml, checkIfExists: true)]), params.format)
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- greenhill
mod('greenhill', 'GREENHILL', 'process_high', 'greenhill:1.1.0--h663a4a6_3',
    """    tuple val(meta), path(hic)
    path(aware)
    path(bubble)
    path(cph)
    path(long)
    val(hic_type)
""",
    """    tuple val(meta), path("${prefix}*"), emit: results
    path "versions.yml"                , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = hic_type == 'paired' ? 'HIC' : 'hic'
    def aware_opt = aware ? "-c ${aware}" + (bubble ? " -b ${bubble.join(' ')}" : '') : ''
    def cph_opt = cph ? "-cph ${cph}" : ''
    def long_opt = long ? "-p ${long.join(' ')}" : ''
    """ + '''
    """
    greenhill \\\\
        ${aware_opt} \\\\
        ${cph_opt} \\\\
        ${long_opt} \\\\
        -t ${task.cpus} \\\\
        ${mode} ${hic} \\\\
        -o ${prefix} \\\\
        ${args}
    """ + versions('greenhill', '1.1.0') + ''',
    """    touch ${meta.id}.scaffolds.fa
    """ + versions('greenhill', '1.1.0'),
    ['greenhill=1.1.0'],
    meta_yml('greenhill', 'Scaffold genomes with GreenHill using Hi-C (and optional long reads / linked reads).',
             ['scaffolding', 'hi-c', 'genome assembly'],
             'greenhill', 'GreenHill: Hi-C scaffolding with long-read support.', 'https://github.com/COMBINE-lab/GreenHill',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'hic': {'type': 'file', 'description': 'Hi-C contact file (pairs or single)'}},
              {'aware': {'type': 'file', 'description': 'Optional linked-read aware file'}},
              {'bubble': {'type': 'file', 'description': 'Optional bubble files accompanying aware'}},
              {'cph': {'type': 'file', 'description': 'Optional CPH contact file (mutually exclusive with aware)'}},
              {'long': {'type': 'file', 'description': 'Optional long reads'}},
              {'hic_type': {'type': 'string', 'description': 'paired (HIC) or single (hic)'}}],
             [{'results': {'type': 'file', 'description': 'GreenHill output files with the given prefix'}}]),
    dict(name='GreenHill', summary='Scaffold a genome with GreenHill using Hi-C contacts, optionally with long or linked reads.',
         params={'hic': P('input', 'file', 'Hi-C contact file.', True),
                 'aware': P('input', 'file', 'Optional linked-read aware file (mutually exclusive with cph).'),
                 'bubble': P('input', 'path', 'Optional directory of bubble files for aware.'),
                 'cph': P('input', 'file', 'Optional CPH contact file.'),
                 'long': P('input', 'path', 'Optional directory of long-read files.'),
                 'hic_type': P('option', 'string', 'Hi-C input type.', enum=['paired', 'single']),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'hic': None, 'aware': None, 'bubble': None, 'cph': None, 'long': None,
                   'hic_type': 'paired', 'args': '', 'outdir': 'results'},
         body="""workflow {
    if (params.aware && params.cph) {
        error "greenhill: 'aware' and 'cph' cannot be used at the same time."
    }
    def files = { p -> p ? file(p, checkIfExists: true).listFiles().findAll { it.isFile() }.toList() : [] }

    GREENHILL(
        Channel.value([[id: 'greenhill'], file(params.hic, checkIfExists: true)]),
        params.aware ? file(params.aware, checkIfExists: true) : [],
        files(params.bubble),
        params.cph ? file(params.cph, checkIfExists: true) : [],
        files(params.long),
        params.hic_type
    )
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- isescan
mod('isescan', 'ISESCAN', 'process_medium', 'isescan:1.7.3--hab16a5f_1',
    """    tuple val(meta), path(contigs)
""",
    """    tuple val(meta), path("*.tsv")  , emit: tsv
    tuple val(meta), path("*.fasta"), emit: fasta, optional: true
    tuple val(meta), path("*.gff")  , emit: gff, optional: true
    path "versions.yml"             , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    """ + '''
    """
    isescan.py \\\\
        --seqfile ${contigs} \\\\
        --output . \\\\
        --nthread ${task.cpus} \\\\
        ${args}
    """ + versions('isescan', '1.7.3') + ''',
    """    touch ${meta.id}.tsv
    """ + versions('isescan', '1.7.3'),
    ['python=3.11', 'isescan=1.7.3'],
    meta_yml('isescan', 'Detect insertion sequence (IS) elements in bacterial genomes with ISEScan.',
             ['insertion sequences', 'transposons', 'bacteria', 'annotation'],
             'isescan', 'ISEScan: IS element detection.', 'https://github.com/xiezhq/ISEScan',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'contigs': {'type': 'file', 'description': 'Genome/contigs FASTA', 'pattern': '*.{fasta,fa,fna}'}}],
             [{'tsv': {'type': 'file', 'description': 'IS element prediction table'}},
              {'fasta': {'type': 'file', 'description': 'IS element sequences'}},
              {'gff': {'type': 'file', 'description': 'IS element GFF annotation'}}]),
    dict(name='ISEScan', summary='Detect insertion sequence (IS) elements in a genome FASTA with ISEScan.',
         params={'contigs': P('input', 'file', 'Genome or contigs FASTA.', True),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'contigs': None, 'args': '', 'outdir': 'results'},
         body="""workflow {
    ISESCAN(Channel.value([[id: 'isescan'], file(params.contigs, checkIfExists: true)]))
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- kneaddata
mod('kneaddata', 'KNEADDATA', 'process_high', 'kneaddata:0.12.4--pyhdfd78af_0',
    """    tuple val(meta), path(reads)
    path(db)
""",
    """    tuple val(meta), path("*_kneaddata.fastq"), emit: reads, optional: true
    tuple val(meta), path("*contam.fastq")    , emit: contaminant, optional: true
    tuple val(meta), path("kneaddata.txt")    , emit: counts
    tuple val(meta), path("*.log")            , emit: log
    path "versions.yml"                       , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def read_list = reads instanceof List ? reads : [reads]
    def input_cmd = read_list.size() == 2 ? "--input1 ${read_list[0]} --input2 ${read_list[1]}" : "--unpaired ${read_list[0]}"
    """ + '''
    """
    kneaddata \\\\
        ${input_cmd} \\\\
        --reference-db ${db} \\\\
        -t ${task.cpus} \\\\
        --output . \\\\
        ${args}

    kneaddata_read_count_table \\\\
        --input . \\\\
        --output kneaddata.txt
    """ + versions('kneaddata', '0.12.4') + ''',
    """    touch sample_kneaddata.fastq
    touch kneaddata.txt
    touch sample.log
    """ + versions('kneaddata', '0.12.4'),
    ['kneaddata=0.12.4'],
    meta_yml('kneaddata', 'Remove host contamination from metagenomic reads with kneaddata.',
             ['metagenomics', 'qc', 'host removal', 'contamination'],
             'kneaddata', 'kneaddata: host-read decontamination.', 'https://github.com/biobakery/kneaddata',
             [{'meta': {'type': 'map', 'description': 'Sample meta map (id, single_end)'}},
              {'reads': {'type': 'file', 'description': 'One (unpaired) or two (paired) FASTQ files'}},
              {'db': {'type': 'directory', 'description': 'kneaddata host reference database directory'}}],
             [{'reads': {'type': 'file', 'description': 'Cleaned FASTQ(s)'}},
              {'contaminant': {'type': 'file', 'description': 'Removed contaminant reads'}},
              {'counts': {'type': 'file', 'description': 'Read count QC table'}},
              {'log': {'type': 'file', 'description': 'kneaddata log'}}]),
    dict(name='kneaddata', summary='Remove host contamination from metagenomic FASTQs with kneaddata and summarise read counts.',
         params={'reads': P('input', 'path_glob', 'FASTQ glob: one file (unpaired) or an R1/R2 pair.', True),
                 'db': P('input', 'path', 'kneaddata host reference database directory.', True),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'reads': None, 'db': None, 'args': '', 'outdir': 'results'},
         body="""workflow {
    reads_ch = Channel.fromFilePairs(params.reads, size: -1)
        .map { sample, reads -> [[id: sample, single_end: reads instanceof List ? reads.size() == 1 : true], reads] }

    KNEADDATA(reads_ch, file(params.db, checkIfExists: true))
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- lefse/format
mod('lefse/format', 'LEFSE_FORMAT', 'process_single', 'lefse:1.1.2--pyhdfd78af_0',
    """    tuple val(meta), path(table)
""",
    """    tuple val(meta), path("${prefix}.lefse"), emit: formatted
    path "versions.yml"                     , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ + '''
    """
    lefse_format_input.py \\\\
        ${table} \\\\
        ${prefix}.lefse \\\\
        ${args}
    """ + versions('lefse', '1.1.2') + ''',
    """    touch ${meta.id}.lefse
    """ + versions('lefse', '1.1.2'),
    ['lefse=1.1.2'],
    meta_yml('lefse_format', 'Format a feature table for LEfSe (lefse_format_input.py).',
             ['metagenomics', 'biomarker', 'differential abundance'],
             'lefse', 'LEfSe: LDA effect size biomarker discovery.', 'https://github.com/SegataLab/lefse',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'table': {'type': 'file', 'description': 'Feature table with class/subclass rows prepared for LEfSe'}}],
             [{'formatted': {'type': 'file', 'description': 'LEfSe-formatted input file'}}]),
    dict(name='lefse format', summary='Format a feature table into LEfSe input (header rows must already be prepared).',
         params={'table': P('input', 'file', 'Feature table with LEfSe class/subclass header rows.', True),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'table': None, 'args': '', 'outdir': 'results'},
         body="""// The calling pipeline prepares the class/subclass/subject header rows
// (the Snakemake wrapper's group_line logic); this module only runs the formatter.
workflow {
    LEFSE_FORMAT(Channel.value([[id: 'lefse'], file(params.table, checkIfExists: true)]))
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- lefse/run
mod('lefse/run', 'LEFSE_RUN', 'process_medium', 'lefse:1.1.2--pyhdfd78af_0',
    """    tuple val(meta), path(formatted)
""",
    """    tuple val(meta), path("${prefix}.res"), emit: results
    path "versions.yml"                   , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ + '''
    """
    lefse_run.py \\\\
        ${formatted} \\\\
        ${prefix}.res \\\\
        ${args}
    """ + versions('lefse', '1.1.2') + ''',
    """    touch ${meta.id}.res
    """ + versions('lefse', '1.1.2'),
    ['lefse=1.1.2'],
    meta_yml('lefse_run', 'Run LEfSe differential abundance analysis (lefse_run.py).',
             ['metagenomics', 'biomarker', 'differential abundance', 'lda'],
             'lefse', 'LEfSe: LDA effect size biomarker discovery.', 'https://github.com/SegataLab/lefse',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'formatted': {'type': 'file', 'description': 'LEfSe-formatted input (from lefse/format)'}}],
             [{'results': {'type': 'file', 'description': 'LEfSe .res results file'}}]),
    dict(name='lefse run', summary='Run LEfSe biomarker discovery on a formatted feature table.',
         params={'formatted': P('input', 'file', 'LEfSe-formatted input file (from lefse format).', True),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'formatted': None, 'args': '', 'outdir': 'results'},
         body="""workflow {
    LEFSE_RUN(Channel.value([[id: 'lefse'], file(params.formatted, checkIfExists: true)]))
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- lefse/res
mod('lefse/res', 'LEFSE_PLOT_RES', 'process_single', 'lefse:1.1.2--pyhdfd78af_0',
    """    tuple val(meta), path(res)
""",
    """    tuple val(meta), path("${prefix}.png"), emit: plot
    path "versions.yml"                   , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ + '''
    """
    lefse_plot_res.py \\\\
        ${res} \\\\
        ${prefix}.png \\\\
        ${args}
    """ + versions('lefse', '1.1.2') + ''',
    """    touch ${meta.id}.png
    """ + versions('lefse', '1.1.2'),
    ['lefse=1.1.2'],
    meta_yml('lefse_plot_res', 'Plot LEfSe LDA effect sizes as a horizontal bar chart.',
             ['metagenomics', 'visualization', 'biomarker'],
             'lefse', 'LEfSe: LDA effect size biomarker discovery.', 'https://github.com/SegataLab/lefse',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'res': {'type': 'file', 'description': 'LEfSe .res results file'}}],
             [{'plot': {'type': 'file', 'description': 'LDA effect size bar plot (PNG)'}}]),
    dict(name='lefse plot_res', summary='Plot LEfSe LDA effect sizes as a bar chart.',
         params={'res': P('input', 'file', 'LEfSe .res results file.', True),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'res': None, 'args': '', 'outdir': 'results'},
         body="""workflow {
    LEFSE_PLOT_RES(Channel.value([[id: 'lefse'], file(params.res, checkIfExists: true)]))
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- lefse/cladogram
mod('lefse/cladogram', 'LEFSE_PLOT_CLADOGRAM', 'process_single', 'lefse:1.1.2--pyhdfd78af_0',
    """    tuple val(meta), path(res)
""",
    """    tuple val(meta), path("${prefix}.png"), emit: plot
    path "versions.yml"                   , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ + '''
    """
    lefse_plot_cladogram.py \\\\
        ${res} \\\\
        ${prefix}.png \\\\
        ${args}
    """ + versions('lefse', '1.1.2') + ''',
    """    touch ${meta.id}.png
    """ + versions('lefse', '1.1.2'),
    ['lefse=1.1.2'],
    meta_yml('lefse_plot_cladogram', 'Plot LEfSe results as a taxonomic cladogram.',
             ['metagenomics', 'visualization', 'biomarker', 'taxonomy'],
             'lefse', 'LEfSe: LDA effect size biomarker discovery.', 'https://github.com/SegataLab/lefse',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'res': {'type': 'file', 'description': 'LEfSe .res results file'}}],
             [{'plot': {'type': 'file', 'description': 'Cladogram figure (PNG)'}}]),
    dict(name='lefse cladogram', summary='Plot LEfSe results as a taxonomic cladogram.',
         params={'res': P('input', 'file', 'LEfSe .res results file.', True),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'res': None, 'args': '', 'outdir': 'results'},
         body="""workflow {
    LEFSE_PLOT_CLADOGRAM(Channel.value([[id: 'lefse'], file(params.res, checkIfExists: true)]))
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- lefse/features
mod('lefse/features', 'LEFSE_PLOT_FEATURES', 'process_single', 'lefse:1.1.2--pyhdfd78af_0',
    """    tuple val(meta), path(dataset), path(res)
""",
    """    tuple val(meta), path("feature_plots/*.png"), emit: plots
    path "versions.yml"                         , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    """ + '''
    """
    mkdir -p feature_plots

    lefse_plot_features.py \\\\
        ${dataset} \\\\
        ${res} \\\\
        feature_plots \\\\
        ${args}
    """ + versions('lefse', '1.1.2') + ''',
    """    mkdir -p feature_plots
    touch feature_plots/feature.png
    """ + versions('lefse', '1.1.2'),
    ['lefse=1.1.2'],
    meta_yml('lefse_plot_features', 'Plot per-feature abundance differences for significant LEfSe biomarkers.',
             ['metagenomics', 'visualization', 'biomarker'],
             'lefse', 'LEfSe: LDA effect size biomarker discovery.', 'https://github.com/SegataLab/lefse',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'dataset': {'type': 'file', 'description': 'LEfSe-formatted input dataset'}},
              {'res': {'type': 'file', 'description': 'LEfSe .res results file'}}],
             [{'plots': {'type': 'file', 'description': 'Per-feature abundance plots (PNG)'}}]),
    dict(name='lefse plot_features', summary='Plot per-feature abundances for significant LEfSe biomarkers.',
         params={'dataset': P('input', 'file', 'LEfSe-formatted input dataset.', True),
                 'res': P('input', 'file', 'LEfSe .res results file.', True),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'dataset': None, 'res': None, 'args': '', 'outdir': 'results'},
         body="""workflow {
    LEFSE_PLOT_FEATURES(Channel.value([
        [id: 'lefse'],
        file(params.dataset, checkIfExists: true),
        file(params.res, checkIfExists: true)
    ]))
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- pysradb
mod('pysradb', 'PYSRADB', 'process_single', 'pysradb:2.5.1--pyhdfd78af_0',
    """    tuple val(meta), val(command), val(query)
""",
    """    tuple val(meta), path("${prefix}*"), emit: results
    path "versions.yml"                , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ + '''
    """
    pysradb \\\\
        ${command} \\\\
        ${query} \\\\
        --saveto ${prefix}.tsv \\\\
        ${args}
    """ + versions('pysradb', '2.5.1') + ''',
    """    touch ${meta.id}.tsv
    """ + versions('pysradb', '2.5.1'),
    ['pysradb=2.5.1'],
    meta_yml('pysradb', 'Query SRA/ENA/DDBJ metadata with pysradb.',
             ['sra', 'metadata', 'download', 'public data'],
             'pysradb', 'pysradb: SRA metadata queries.', 'https://github.com/saketkc/pysradb',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'command': {'type': 'string', 'description': 'pysradb subcommand (e.g. metadata, search, download)'}},
              {'query': {'type': 'string', 'description': 'Accession or search query'}}],
             [{'results': {'type': 'file', 'description': 'Saved query result table'}}]),
    dict(name='pysradb', summary='Query SRA metadata (or search/download) with pysradb and save the result table.',
         params={'command': P('option', 'string', 'pysradb subcommand, e.g. metadata, search.', True),
                 'query': P('option', 'string', 'Accession or query string, e.g. SRP000001.', True),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'command': 'metadata', 'query': None, 'args': '', 'outdir': 'results'},
         body="""// Requires network access to NCBI/ENA.
workflow {
    PYSRADB(Channel.value([[id: 'pysradb'], params.command, params.query]))
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))

# ---------------------------------------------------------------- virsorter2
mod('virsorter2', 'VIRSORTER2', 'process_high', 'virsorter:2.2.4--pyhdfd78af_2',
    """    tuple val(meta), path(contigs)
    path(db_dir)
    val(mode)
    val(include_groups)
    val(min_length)
""",
    """    tuple val(meta), path("results/final-viral-combined.fa"), emit: combined, optional: true
    tuple val(meta), path("results/final-viral-score.tsv")  , emit: score, optional: true
    tuple val(meta), path("results/for-dramv")              , emit: dramv, optional: true
    path "versions.yml"                                   , emit: versions, topic: versions
""",
    """    def args = task.ext.args ?: ''
    def groups_opt = include_groups ? "--include-groups ${include_groups}" : ''
    def db_opt = db_dir ? "--db-dir ${db_dir}" : ''
    def setup = ''
    if (db_dir) {
        setup = "if [ ! -f ${db_dir}/.vs2_setup.done ]; then virsorter setup -d ${db_dir} -j ${task.cpus} && touch ${db_dir}/.vs2_setup.done; fi; virsorter config --init-source --db-dir=${db_dir}"
    }
    """ + '''
    """
    ${setup}

    virsorter run \\\\
        -w results \\\\
        -i ${contigs} \\\\
        --min-length ${min_length} \\\\
        -j ${task.cpus} \\\\
        ${groups_opt} \\\\
        ${db_opt} \\\\
        ${args} \\\\
        ${mode}
    """ + versions('virsorter', '2.2.4') + ''',
    """    mkdir -p results/for-dramv
    touch results/final-viral-combined.fa
    touch results/final-viral-score.tsv
    """ + versions('virsorter', '2.2.4'),
    ['virsorter=2.2.4'],
    meta_yml('virsorter2', 'Identify viral sequences from metagenomic assemblies with VirSorter2.',
             ['virus', 'metagenomics', 'classification'],
             'virsorter', 'VirSorter2: viral sequence detection.', 'https://github.com/jiarong/VirSorter2',
             [{'meta': {'type': 'map', 'description': 'Sample meta map'}},
              {'contigs': {'type': 'file', 'description': 'Assembly contigs FASTA'}},
              {'db_dir': {'type': 'directory', 'description': 'VirSorter2 database directory (set up on first use)'}},
              {'mode': {'type': 'string', 'description': 'all (detect + classify) or classify'}},
              {'include_groups': {'type': 'string', 'description': 'Comma-separated viral groups to include'}},
              {'min_length': {'type': 'integer', 'description': 'Minimum contig length'}}],
             [{'combined': {'type': 'file', 'description': 'Final viral combined FASTA'}},
              {'score': {'type': 'file', 'description': 'Final viral score table'}},
              {'dramv': {'type': 'directory', 'description': 'for-dramv inputs for DRAM-v'}}]),
    dict(name='VirSorter2', summary='Detect and classify viral contigs with VirSorter2 (produces for-dramv inputs for DRAM-v).',
         params={'contigs': P('input', 'file', 'Assembly contigs FASTA.', True),
                 'db_dir': P('input', 'path', 'VirSorter2 database directory; set up automatically on first use.'),
                 'mode': P('option', 'string', 'Run mode.', enum=['all', 'classify']),
                 'include_groups': P('option', 'string', 'Comma-separated viral groups to include (e.g. dsDNAphage,ssDNA).'),
                 'min_length': P('option', 'integer', 'Minimum contig length.'),
                 'args': ARGS, 'outdir': OUTDIR},
         defaults={'contigs': None, 'db_dir': None, 'mode': 'all', 'include_groups': None,
                   'min_length': 1500, 'args': '', 'outdir': 'results'},
         body="""workflow {
    VIRSORTER2(
        Channel.value([[id: 'virsorter2'], file(params.contigs, checkIfExists: true)]),
        params.db_dir ? file(params.db_dir, checkIfExists: true) : [],
        params.mode,
        params.include_groups ?: '',
        params.min_length
    )
}
""",
         extargs="        ext.args = { params.args ?: '' }\n"))


# ================================================================ emit files

PROC_HEADER = """process {proc} {{
    tag "$meta.id"
    label '{label}'

    conda "${{moduleDir}}/environment.yml"
{container}
    input:
{inputs}
    output:
{outputs}
    when:
    task.ext.when == null || task.ext.when

    script:
{script}
    stub:
{stub}
}}
"""

WRAPPER_HEADER = """#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/{src}, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include {{ {proc} }} from '../main.nf'

"""

CONFIG = """// Launcher config for the {proc} wrapper adapter. Not part of the fixed
// wrapper/{{main.nf,params.json,wrapper.yaml}} contract — Nextflow picks this
// up automatically because it sits beside wrapper/main.nf.

process {{
    cpus   = {cpus}
    memory = '{mem}'
    time   = '{time}'

    withName: '{proc}' {{
        publishDir = [
            path: {{ "${{params.outdir}}/{outdir_name}" }},
            mode: 'copy'
        ]
{extargs}    }}
}}

profiles {{
    docker {{
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
        conda.enabled = false
        singularity.enabled = false
    }}
    singularity {{
        singularity.enabled = true
        singularity.autoMounts = true
        docker.enabled = false
        conda.enabled = false
    }}
    conda {{
        conda.enabled = true
        docker.enabled = false
        singularity.enabled = false
    }}
}}
"""

def main():
    assert len(MODULES) == 15, len(MODULES)
    for rel, m in MODULES.items():
        base = os.path.join(DST, rel)
        wdir = os.path.join(base, 'wrapper')
        if os.path.exists(os.path.join(wdir, 'wrapper.yaml')):
            print('skip (exists)', rel)
            continue
        os.makedirs(wdir, exist_ok=True)
        # module main.nf
        has_meta = 'val(meta)' in m['inputs']
        script = re.sub(r'"""\s*\+\s*versions\(\'([^\']+)\',\s*\'([^\']+)\'\)\s*\+\s*$',
                        lambda mm: versions(mm.group(1), mm.group(2)).rstrip() + '\n    """', m['script'])
        stub = '    """\n' + m['stub'] + '\n    """'
        header = PROC_HEADER.format(proc=m['proc'], label=m['label'], container=container(m['tag']),
                                    inputs=m['inputs'], outputs=m['outputs'], script=script, stub=stub)
        if not has_meta:
            header = header.replace('    tag "$meta.id"\n', '')
        open(os.path.join(base, 'main.nf'), 'w').write(header)
        # meta.yml
        with open(os.path.join(base, 'meta.yml'), 'w') as f:
            yaml.safe_dump(m['meta'], f, sort_keys=False, allow_unicode=True, width=100)
        # environment.yml
        env = {'channels': ['conda-forge', 'bioconda', 'nodefaults'],
               'dependencies': m['env_deps']}
        with open(os.path.join(base, 'environment.yml'), 'w') as f:
            yaml.safe_dump(env, f, sort_keys=False)
        # wrapper
        w = m['wrapper']
        decl = ''.join(f'params.{k} = null\n' for k in w['params'])
        src = rel
        open(os.path.join(wdir, 'main.nf'), 'w').write(
            WRAPPER_HEADER.format(proc=m['proc'], src=src) + decl + '\n' + w['body'])
        y = {'id': 'local/modules/' + rel.replace('/', '-').replace('_', '-'),
             'name': w['name'], 'summary': w['summary'], 'params': w['params'],
             'outputs': {'output': {'type': 'directory', 'path': '${outdir}/' + m['proc'].lower(), 'primary': True}}}
        with open(os.path.join(wdir, 'wrapper.yaml'), 'w') as f:
            f.write('# Agent-facing wrapper for the local ' + rel + ' module at ../main.nf,\n')
            f.write('# ported from the Snakemake Workflow project (no nf-core equivalent).\n')
            f.write('# Follows docs/design/phi-wrapper-agent-composition-design.md section 3.\n')
            yaml.safe_dump(y, f, sort_keys=False, allow_unicode=True, width=100)
        assert set(w['defaults']) == set(w['params']), rel
        open(os.path.join(wdir, 'params.json'), 'w').write(json.dumps(w['defaults'], indent=2) + '\n')
        cpus, mem, time = RESOURCES[m['label']]
        open(os.path.join(wdir, 'nextflow.config'), 'w').write(
            CONFIG.format(proc=m['proc'], cpus=cpus, mem=mem, time=time,
                          outdir_name=m['proc'].lower(), extargs=w['extargs']))
        print('ported', rel)

if __name__ == '__main__':
    main()