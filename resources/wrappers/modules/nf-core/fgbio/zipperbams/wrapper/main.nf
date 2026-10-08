#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { FGBIO_ZIPPERBAMS } from '../main.nf'


params.mapped     = null
params.unmapped   = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', single_end: false], file(params.mapped, checkIfExists: true), file(params.unmapped, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true), file(params.dict, checkIfExists: true)])

    FGBIO_ZIPPERBAMS(in_ch, fasta_ch)
}
