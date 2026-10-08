#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_APPLYVQSR } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.recal      = null
params.recal_idx  = null
params.tranches   = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'applied'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true),
                            file(params.recal, checkIfExists: true), file(params.recal_idx, checkIfExists: true), file(params.tranches, checkIfExists: true)])
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))

    GATK4_APPLYVQSR(vcf_ch, fasta_ch, fai_ch, dict_ch)
}
