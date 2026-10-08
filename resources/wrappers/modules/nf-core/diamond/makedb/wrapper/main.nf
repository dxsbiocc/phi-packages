#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../main.nf'

params.fasta  = null
params.outdir = null

workflow {
    fasta_ch = Channel.value([[id: 'diamonddb'], file(params.fasta, checkIfExists: true)])
    empty_ch = Channel.value([])

    DIAMOND_MAKEDB(fasta_ch, empty_ch, empty_ch, empty_ch)
}
