#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_CRAMSIZE } from '../main.nf'

params.cram   = null
params.outdir = null

workflow {
    cram_ch = Channel.value([[id: 'test'], file(params.cram, checkIfExists: true)])

    SAMTOOLS_CRAMSIZE(cram_ch)
}
