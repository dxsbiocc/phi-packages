#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SEQKIT_SLIDING } from '../main.nf'


params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], file(params.fasta, checkIfExists: true)])

    SEQKIT_SLIDING(in_ch)
}
