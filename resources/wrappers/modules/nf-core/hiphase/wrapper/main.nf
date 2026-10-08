#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { HIPHASE } from '../main.nf'


params.bam        = null
params.bai        = null
params.vcf        = null
params.csi        = null
params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], [file(params.bam, checkIfExists: true)], [file(params.bai, checkIfExists: true)],
                          [file(params.vcf, checkIfExists: true)], [file(params.csi, checkIfExists: true)], [], [], []])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true), []])

    HIPHASE(in_ch, fasta_ch, false, false, false, false, false, 'tsv')
}
