#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_VARIANTSTOTABLE } from '../main.nf'


params.vcf        = null
params.idx        = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'test'], file(params.vcf, checkIfExists: true), file(params.idx, checkIfExists: true), [], [], []])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])

    GATK4_VARIANTSTOTABLE(vcf_ch, fasta_ch, fai_ch, dict_ch)
}
