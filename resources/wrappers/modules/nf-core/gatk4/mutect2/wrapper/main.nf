#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_MUTECT2 } from '../main.nf'


params.tumour_bam = null
params.tumour_bai = null
params.normal_bam = null
params.normal_bai = null
params.germline   = null
params.germline_tbi = null
params.pon        = null
params.pon_tbi    = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    input_ch = Channel.value([[id: 'test', normal_id: 'normal', tumor_id: 'tumour'],
                              [file(params.normal_bam, checkIfExists: true), file(params.tumour_bam, checkIfExists: true)],
                              [file(params.normal_bai, checkIfExists: true), file(params.tumour_bai, checkIfExists: true)], []])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], [file(params.fai, checkIfExists: true)], []])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])
    empty_ch = Channel.value([])

    GATK4_MUTECT2(input_ch, fasta_ch, fai_ch, dict_ch, empty_ch, empty_ch,
                  Channel.value(file(params.germline, checkIfExists: true)), Channel.value(file(params.germline_tbi, checkIfExists: true)),
                  Channel.value(file(params.pon, checkIfExists: true)), Channel.value(file(params.pon_tbi, checkIfExists: true)))
}
