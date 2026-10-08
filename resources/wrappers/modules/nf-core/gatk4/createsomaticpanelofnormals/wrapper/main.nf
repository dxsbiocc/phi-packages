#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_CREATESOMATICPANELOFNORMALS } from '../main.nf'
include { UNTAR } from '../../../untar/main.nf'

params.genomicsdb = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    db_ch = Channel.value([[id: 'test'], file(params.genomicsdb, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])
    UNTAR(db_ch)
    GATK4_CREATESOMATICPANELOFNORMALS(UNTAR.out.untar, fasta_ch, fai_ch, dict_ch)
}
