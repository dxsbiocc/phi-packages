#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_SORTVCF } from '../main.nf'


params.vcf        = null
params.fasta      = null
params.dict       = null
params.outdir     = null

workflow {
    vcf_ch   = Channel.value([[id: 'test'], file(params.vcf, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])

    PICARD_SORTVCF(vcf_ch, fasta_ch, dict_ch)
}
