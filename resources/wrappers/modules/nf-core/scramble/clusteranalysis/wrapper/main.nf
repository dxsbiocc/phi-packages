#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SCRAMBLE_CLUSTERANALYSIS } from '../main.nf'
include { SCRAMBLE_CLUSTERIDENTIFIER } from '../../clusteridentifier/main.nf'

params.bam        = null
params.bai        = null
params.fasta      = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)])

    SCRAMBLE_CLUSTERIDENTIFIER(bam_ch, [[], []])
    SCRAMBLE_CLUSTERANALYSIS(SCRAMBLE_CLUSTERIDENTIFIER.out.clusters, fasta_ch, [])
}
