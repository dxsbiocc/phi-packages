#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_BASERECALIBRATOR } from '../main.nf'

params.bam         = null
params.bai         = null
params.fasta       = null
params.fai         = null
params.dict        = null
params.known_sites = null
params.known_sites_tbi = null
params.outdir      = null

workflow {
    input_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), []])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])
    known_ch = Channel.value([[id: 'known'], file(params.known_sites, checkIfExists: true)])
    known_tbi_ch = Channel.value([[id: 'known'], file(params.known_sites_tbi, checkIfExists: true)])

    GATK4_BASERECALIBRATOR(input_ch, fasta_ch, fai_ch, dict_ch, known_ch, known_tbi_ch)
}
