#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GRIDSS_GRIDSS } from '../main.nf'
include { BWA_INDEX } from '../../../bwa/index/main.nf'

params.bam        = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'fasta_fai'], file(params.fai, checkIfExists: true)])
    BWA_INDEX(fasta_ch)

    GRIDSS_GRIDSS(in_ch, fasta_ch, fai_ch, BWA_INDEX.out.index)
}
