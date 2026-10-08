#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_COLLECTSVEVIDENCE } from '../main.nf'


params.cram       = null
params.crai       = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    cram_ch = Channel.value([[id: 'test', single_end: false], file(params.cram, checkIfExists: true), file(params.crai, checkIfExists: true), [], []])
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))

    GATK4_COLLECTSVEVIDENCE(cram_ch, fasta_ch, fai_ch, dict_ch)
}
