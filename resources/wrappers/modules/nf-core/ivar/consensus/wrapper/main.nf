#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { IVAR_CONSENSUS } from '../main.nf'


params.bam        = null
params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true)])

    IVAR_CONSENSUS(in_ch, file(params.fasta, checkIfExists: true), false)
}
