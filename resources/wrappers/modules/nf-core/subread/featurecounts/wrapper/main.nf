#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
//
// `merge: true` counts every BAM in one featureCounts call — featureCounts
// itself then writes one column per BAM — and reshapes that table into a
// plain gene x sample matrix (`gene_id` + one column per sample, named after
// the BAM's simpleName) that differential-expression wrappers read directly.
nextflow.enable.dsl = 2

include { SUBREAD_FEATURECOUNTS } from '../main.nf'

params.bam          = null
params.gtf          = null
params.outdir       = null
params.strandedness = 'unstranded'
params.feature_type = 'exon'
params.single_end   = false
params.merge        = false

// featureCounts' own table: a `#` command line, then Geneid, Chr, Start, End,
// Strand, Length and one count column per input BAM (named by its path).
process FEATURECOUNTS_MATRIX {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(table)

    output:
    tuple val(meta), path("${meta.id}.counts_matrix.tsv"), emit: matrix

    script:
    """
    awk 'BEGIN { FS = OFS = "\\t" }
         /^#/ { next }
         !header {
             printf "gene_id"
             for (i = 7; i <= NF; i++) {
                 name = \$i
                 sub(/.*\\//, "", name)
                 sub(/\\..*/, "", name)
                 printf "%s%s", OFS, name
             }
             print ""
             header = 1
             next
         }
         {
             printf "%s", \$1
             for (i = 7; i <= NF; i++) printf "%s%s", OFS, \$i
             print ""
         }' ${table} > ${meta.id}.counts_matrix.tsv
    """
}

def metaFor(id) {
    [id: id, single_end: params.single_end, strandedness: params.strandedness]
}

workflow {
    gtf = file(params.gtf, checkIfExists: true)
    bams = Channel.fromPath(params.bam, checkIfExists: true)

    if (params.merge) {
        counts_ch = bams
            .toSortedList { a, b -> a.name <=> b.name }
            .map { files -> [metaFor('merged'), files, gtf] }
    } else {
        counts_ch = bams.map { bam -> [metaFor(bam.simpleName), bam, gtf] }
    }

    SUBREAD_FEATURECOUNTS(counts_ch)

    if (params.merge) {
        FEATURECOUNTS_MATRIX(SUBREAD_FEATURECOUNTS.out.counts)
    }
}
