#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMBAMBA_DEPTH } from '../main.nf'


params.bam        = null
params.bai        = null
params.mode       = 'base'
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'sample', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])

    SAMBAMBA_DEPTH(bam_ch, [[], []], params.mode)
}
