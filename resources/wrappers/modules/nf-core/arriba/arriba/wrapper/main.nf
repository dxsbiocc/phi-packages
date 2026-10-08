#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { ARRIBA_ARRIBA } from '../main.nf'
include { STAR_GENOMEGENERATE } from '../../../star/genomegenerate/main.nf'
include { STAR_ALIGN } from '../../../star/align/main.nf'

params.reads1     = null
params.reads2     = null
params.fasta      = null
params.gtf        = null
params.outdir     = null

workflow {
    fasta_ch = Channel.value([[id: 'fasta'], [file(params.fasta, checkIfExists: true)]])
    gtf_ch   = Channel.value([[id: 'gtf'], [file(params.gtf, checkIfExists: true)]])
    STAR_GENOMEGENERATE(fasta_ch, gtf_ch)
    reads_ch = Channel.value([[id: 'sample', single_end: false], [file(params.reads1, checkIfExists: true), file(params.reads2, checkIfExists: true)]])
    STAR_ALIGN(reads_ch, STAR_GENOMEGENERATE.out.index, gtf_ch, false)

    ARRIBA_ARRIBA(STAR_ALIGN.out.bam, fasta_ch, gtf_ch, [], [], [], [])
}
