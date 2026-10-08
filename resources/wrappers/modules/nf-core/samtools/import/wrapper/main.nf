#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_IMPORT } from '../main.nf'

params.reads  = null
params.outdir = null

workflow {
    reads_ch = Channel.value([[id: 'test', single_end: true], file(params.reads, checkIfExists: true)])

    SAMTOOLS_IMPORT(reads_ch)
}
