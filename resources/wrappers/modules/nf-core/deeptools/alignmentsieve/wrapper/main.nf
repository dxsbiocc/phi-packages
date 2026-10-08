#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { DEEPTOOLS_ALIGNMENTSIEVE } from '../main.nf'

params.bam    = null
params.bai    = null
params.outdir = null

workflow {
    bam_ch = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])

    DEEPTOOLS_ALIGNMENTSIEVE(bam_ch)
}
