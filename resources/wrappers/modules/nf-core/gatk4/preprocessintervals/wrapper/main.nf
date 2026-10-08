#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_PREPROCESSINTERVALS } from '../main.nf'


params.fasta      = null
params.fai        = null
params.dict       = null
params.exclude    = null
params.outdir     = null

workflow {
    fasta_ch   = Channel.value([[id: 'test'], file(params.fasta, checkIfExists: true)])
    fai_ch     = Channel.value([[id: 'test'], file(params.fai, checkIfExists: true)])
    dict_ch    = Channel.value([[id: 'test'], file(params.dict, checkIfExists: true)])
    include_ch = Channel.value([[], []])
    exclude_ch = Channel.value([[id: 'test'], file(params.exclude, checkIfExists: true)])

    GATK4_PREPROCESSINTERVALS(fasta_ch, fai_ch, dict_ch, include_ch, exclude_ch)
}
