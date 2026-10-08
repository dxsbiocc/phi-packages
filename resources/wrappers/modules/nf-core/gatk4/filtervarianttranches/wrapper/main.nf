#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_FILTERVARIANTTRANCHES } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.resource   = null
params.resource_tbi = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'test'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true), []])
    res_ch = Channel.value(file(params.resource, checkIfExists: true))
    resi_ch = Channel.value(file(params.resource_tbi, checkIfExists: true))
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))

    GATK4_FILTERVARIANTTRANCHES(vcf_ch, res_ch, resi_ch, fasta_ch, fai_ch, dict_ch)
}
