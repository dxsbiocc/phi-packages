#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_COLLATEFASTQ } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    bam_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true)])
    fasta_ch = Channel.value([[:], [], []])

    SAMTOOLS_COLLATEFASTQ(bam_ch, fasta_ch, false)
}
