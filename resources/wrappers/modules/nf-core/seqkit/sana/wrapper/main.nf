#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SEQKIT_SANA } from '../main.nf'


params.reads      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result', single_end: true], file(params.reads, checkIfExists: true)])

    SEQKIT_SANA(in_ch)
}
