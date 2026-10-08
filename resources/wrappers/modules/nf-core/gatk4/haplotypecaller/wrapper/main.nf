#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_HAPLOTYPECALLER } from '../main.nf'

params.bam    = null
params.bai    = null
params.fasta  = null
params.fai    = null
params.dict   = null
params.outdir = null

workflow {
    input_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), [], []])

    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])
    empty_ch = Channel.value([[:], []])

    GATK4_HAPLOTYPECALLER(input_ch, fasta_ch, fai_ch, dict_ch, empty_ch, empty_ch)
}
