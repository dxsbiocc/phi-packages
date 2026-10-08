#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// Composes the vendored ucsc/gtftogenepred module (../../../ucsc/gtftogenepred/main.nf)
// to build the refFlat annotation Picard needs from a GTF.
nextflow.enable.dsl = 2

include { PICARD_COLLECTRNASEQMETRICS } from '../main.nf'
include { UCSC_GTFTOGENEPRED } from '../../../ucsc/gtftogenepred/main.nf'

params.bam        = null
params.gtf        = null
params.fasta      = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'test', single_end: false, strandedness: 'forward'], file(params.bam, checkIfExists: true)])
    gtf_ch = Channel.value([[id: 'genome'], [file(params.gtf, checkIfExists: true)]])

    UCSC_GTFTOGENEPRED(gtf_ch)

    refflat_ch = UCSC_GTFTOGENEPRED.out.refflat.map { meta, f -> f }
    fasta_ch   = Channel.value(file(params.fasta, checkIfExists: true))

    PICARD_COLLECTRNASEQMETRICS(bam_ch, refflat_ch, fasta_ch, Channel.value([]))
}
