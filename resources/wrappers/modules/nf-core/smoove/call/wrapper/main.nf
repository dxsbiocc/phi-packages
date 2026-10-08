#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SMOOVE_CALL } from '../main.nf'


params.bam        = null
params.bai        = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), []])
    fasta_ch = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'fai'], file(params.fai, checkIfExists: true)])

    SMOOVE_CALL(in_ch, fasta_ch, fai_ch)
}
