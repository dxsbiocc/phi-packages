#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SYLPH_SKETCH } from '../main.nf'

params.reads     = null
params.reference = null
params.outdir    = null

workflow {
    reads_ch = Channel.value([[id: 'test', single_end: true], [file(params.reads, checkIfExists: true)]])
    reference_ch = Channel.value(file(params.reference, checkIfExists: true))

    SYLPH_SKETCH(reads_ch, reference_ch)
}
