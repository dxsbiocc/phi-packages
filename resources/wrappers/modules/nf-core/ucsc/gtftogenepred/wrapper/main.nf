#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { UCSC_GTFTOGENEPRED } from '../main.nf'

params.gtf    = null
params.outdir = null

workflow {
    gtf_ch = Channel.value([[id: 'test'], [file(params.gtf, checkIfExists: true)]])

    UCSC_GTFTOGENEPRED(gtf_ch)
}
