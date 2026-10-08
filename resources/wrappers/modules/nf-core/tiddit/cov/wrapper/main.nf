#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { TIDDIT_COV } from '../main.nf'


params.cram       = null
params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', single_end: false], file(params.cram, checkIfExists: true), []])
    fasta_ch = Channel.value([[:], file(params.fasta, checkIfExists: true)])

    TIDDIT_COV(in_ch, fasta_ch)
}
