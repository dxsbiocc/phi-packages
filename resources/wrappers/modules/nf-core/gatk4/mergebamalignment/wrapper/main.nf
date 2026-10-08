#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_MERGEBAMALIGNMENT } from '../main.nf'


params.aligned    = null
params.unmapped   = null
params.fasta      = null
params.dict       = null
params.outdir     = null

workflow {
    bam_ch   = Channel.value([[id: 'merged'], file(params.aligned, checkIfExists: true), file(params.unmapped, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])

    GATK4_MERGEBAMALIGNMENT(bam_ch, fasta_ch, dict_ch)
}
