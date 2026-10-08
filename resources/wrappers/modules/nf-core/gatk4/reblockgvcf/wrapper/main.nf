#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_REBLOCKGVCF } from '../main.nf'


params.gvcf       = null
params.tbi        = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    gvcf_ch  = Channel.value([[id: 'test', single_end: false], file(params.gvcf, checkIfExists: true), file(params.tbi, checkIfExists: true), []])
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))
    empty_ch = Channel.value([])

    GATK4_REBLOCKGVCF(gvcf_ch, fasta_ch, fai_ch, dict_ch, empty_ch, empty_ch)
}
