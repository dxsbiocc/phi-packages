#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_APPLYBQSR } from '../main.nf'

params.bam    = null
params.bai    = null
params.table  = null
params.fasta  = null
params.fai    = null
params.dict   = null
params.outdir = null

workflow {
    input_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), file(params.table, checkIfExists: true), []])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true), file(params.dict, checkIfExists: true)])

    GATK4_APPLYBQSR(input_ch, fasta_ch, 'bam')
}
