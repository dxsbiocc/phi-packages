#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEDTOOLS_GETFASTA } from '../main.nf'


params.bed        = null
params.fasta      = null
params.outdir     = null

workflow {
    bed_ch = Channel.value([[id: 'result'], file(params.bed, checkIfExists: true)])

    BEDTOOLS_GETFASTA(bed_ch, file(params.fasta, checkIfExists: true))
}
