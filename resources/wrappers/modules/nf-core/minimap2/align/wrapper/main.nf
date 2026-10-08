#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { MINIMAP2_ALIGN } from '../main.nf'


params.reads      = null
params.fasta      = null
params.outdir     = null

workflow {
    reads_ch = Channel.value([[id: 'sample', single_end: true], file(params.reads, checkIfExists: true)])
    ref_ch   = Channel.value([[id: 'ref'], file(params.fasta, checkIfExists: true)])

    MINIMAP2_ALIGN(reads_ch, ref_ch, true, [], false, false)
}
