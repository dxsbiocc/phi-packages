#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PBSV_DISCOVER } from '../main.nf'


params.bam        = null
params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'ref'], file(params.fasta, checkIfExists: true)])

    PBSV_DISCOVER(in_ch, fasta_ch)
}
