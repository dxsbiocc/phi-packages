#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_ADDREPLACERG } from '../main.nf'

params.bam        = null
params.read_group = "'@RG\\tID:1\\tLB:lib1\\tPL:ILLUMINA\\tSM:test\\tPU:barcode1'"
params.outdir     = null

workflow {
    bam_ch   = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true), [], params.read_group])
    fasta_ch = Channel.value([[:], [], [], []])

    SAMTOOLS_ADDREPLACERG(bam_ch, fasta_ch)
}
