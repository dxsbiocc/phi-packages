#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { TRUVARI_BENCH } from '../main.nf'


params.sv         = null
params.sv_tbi     = null
params.sv2        = null
params.sv2_tbi    = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.sv, checkIfExists: true), file(params.sv_tbi, checkIfExists: true),
                           file(params.sv2, checkIfExists: true), file(params.sv2_tbi, checkIfExists: true), []])
    fasta_ch = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'fai'], file(params.fai, checkIfExists: true)])

    TRUVARI_BENCH(in_ch, fasta_ch, fai_ch)
}
