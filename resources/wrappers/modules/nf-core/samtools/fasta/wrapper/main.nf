#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_FASTA } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    bam_ch = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true)])

    SAMTOOLS_FASTA(bam_ch, false)
}
