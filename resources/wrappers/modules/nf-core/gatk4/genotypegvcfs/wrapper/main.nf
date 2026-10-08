#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_GENOTYPEGVCFS } from '../main.nf'

params.gvcf       = null
params.gvcf_index = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    gvcf_ch  = Channel.value([[id: 'test'], file(params.gvcf, checkIfExists: true), file(params.gvcf_index, checkIfExists: true), [], []])
    fasta_ch = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'fasta'], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'fasta'], file(params.dict, checkIfExists: true)])
    empty_ch = Channel.value([[:], []])

    GATK4_GENOTYPEGVCFS(gvcf_ch, fasta_ch, fai_ch, dict_ch, empty_ch, empty_ch)
}
