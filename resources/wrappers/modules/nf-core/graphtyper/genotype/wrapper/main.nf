#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GRAPHTYPER_GENOTYPE } from '../main.nf'


params.bam        = null
params.bai        = null
params.fasta      = null
params.fai        = null
params.regions    = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'ref'], file(params.fasta, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'ref_index'], file(params.fai, checkIfExists: true)])

    GRAPHTYPER_GENOTYPE(bam_ch, fasta_ch, fai_ch, file(params.regions, checkIfExists: true))
}
