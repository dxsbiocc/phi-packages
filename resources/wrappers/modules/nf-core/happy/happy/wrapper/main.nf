#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { HAPPY_HAPPY } from '../main.nf'


params.query      = null
params.truth      = null
params.regions_bed = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.query, checkIfExists: true), file(params.truth, checkIfExists: true), file(params.regions_bed, checkIfExists: true), []])
    fasta_ch = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'fai'], file(params.fai, checkIfExists: true)])
    empty_ch = Channel.value([[], []])

    HAPPY_HAPPY(in_ch, fasta_ch, fai_ch, empty_ch, empty_ch, empty_ch)
}
