#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BBMAP_BBNORM } from '../main.nf'


params.reads      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result', single_end: true], file(params.reads, checkIfExists: true)])

    BBMAP_BBNORM(in_ch)
}
