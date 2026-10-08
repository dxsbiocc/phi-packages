#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GRIDSS_PREPROCESS } from '../main.nf'
include { BWA_INDEX } from '../../../bwa/index/main.nf'

params.bam        = null
params.bai        = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    bam_ch = Channel.of([[id: 'sample'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    fasta_ch = Channel.of([[id: 'fasta'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])
    BWA_INDEX(Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)]))
    ref_ch = fasta_ch.join(BWA_INDEX.out.index)

    GRIDSS_PREPROCESS(bam_ch, ref_ch)
}
