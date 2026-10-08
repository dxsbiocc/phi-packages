#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_CONSENSUS } from '../main.nf'

params.bam    = null
params.bai    = null
params.outdir = null

workflow {
    bam_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])

    SAMTOOLS_CONSENSUS(bam_ch)
}
