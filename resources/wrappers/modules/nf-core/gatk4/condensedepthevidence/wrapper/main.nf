#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_CONDENSEDEPTHEVIDENCE } from '../main.nf'


params.depth      = null
params.depth_idx  = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    depth_ch = Channel.value([[id: 'test', single_end: false], file(params.depth, checkIfExists: true), file(params.depth_idx, checkIfExists: true)])
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))

    GATK4_CONDENSEDEPTHEVIDENCE(depth_ch, fasta_ch, fai_ch, dict_ch)
}
