#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BBMAP_PILEUP } from '../main.nf'


params.bam        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result', single_end: false], file(params.bam, checkIfExists: true)])

    BBMAP_PILEUP(in_ch)
}
