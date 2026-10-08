#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { LONGPHASE_PHASE } from '../main.nf'


params.bam        = null
params.bai        = null
params.vcf        = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), file(params.vcf, checkIfExists: true), [], []])
    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'reference'], file(params.fai, checkIfExists: true)])

    LONGPHASE_PHASE(in_ch, fasta_ch, fai_ch)
}
