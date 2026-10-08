#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { MANTA_SOMATIC } from '../main.nf'


params.normal     = null
params.normal_idx = null
params.tumor      = null
params.tumor_idx  = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', single_end: false], file(params.normal, checkIfExists: true), file(params.normal_idx, checkIfExists: true),
                           file(params.tumor, checkIfExists: true), file(params.tumor_idx, checkIfExists: true), [], []])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])
    config_ch = Channel.of('[manta]', 'enableRemoteReadRetrievalForInsertionsInGermlineCallingModes = 0')
                       .collectFile(name: 'manta_options.ini', newLine: true)

    MANTA_SOMATIC(in_ch, fasta_ch, fai_ch, config_ch)
}
