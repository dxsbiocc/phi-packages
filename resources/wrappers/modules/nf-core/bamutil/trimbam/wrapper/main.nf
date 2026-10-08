#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BAMUTIL_TRIMBAM } from '../main.nf'


params.bam        = null
params.trim_left  = 2
params.trim_right = 2
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'sample', single_end: false], file(params.bam, checkIfExists: true), params.trim_left, params.trim_right])

    BAMUTIL_TRIMBAM(bam_ch)
}
