#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { MMSEQS_EASYCLUSTER } from '../main.nf'


params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result', single_end: false], file(params.fasta, checkIfExists: true)])

    MMSEQS_EASYCLUSTER(in_ch)
}
