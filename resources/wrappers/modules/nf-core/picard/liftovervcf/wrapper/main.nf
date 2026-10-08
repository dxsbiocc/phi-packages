#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_LIFTOVERVCF } from '../main.nf'


params.vcf        = null
params.dict       = null
params.fasta      = null
params.chain      = null
params.outdir     = null

workflow {
    vcf_ch   = Channel.value([[id: 'test'], file(params.vcf, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    chain_ch = Channel.value([[id: 'genome'], file(params.chain, checkIfExists: true)])

    PICARD_LIFTOVERVCF(vcf_ch, dict_ch, fasta_ch, chain_ch)
}
