#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_EXTRACTFINGERPRINT } from '../main.nf'


params.bam        = null
params.bai        = null
params.haplotype_map = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])

    PICARD_EXTRACTFINGERPRINT(bam_ch, Channel.value(file(params.haplotype_map, checkIfExists: true)), Channel.value(file(params.fasta, checkIfExists: true)),
                              Channel.value(file(params.fai, checkIfExists: true)), Channel.value(file(params.dict, checkIfExists: true)))
}
