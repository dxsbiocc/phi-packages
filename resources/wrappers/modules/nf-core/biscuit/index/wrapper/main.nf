#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BISCUIT_INDEX } from '../main.nf'


params.fasta      = null
params.outdir     = null

workflow {
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])

    BISCUIT_INDEX(fasta_ch)
}
