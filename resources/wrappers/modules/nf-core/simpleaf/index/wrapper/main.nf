#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Only the genome fasta + gtf input path is exposed (the module's other
// three input slots — transcript fasta, probe csv, feature csv — are
// mutually exclusive alternatives; this is the common case).
nextflow.enable.dsl = 2

include { SIMPLEAF_INDEX } from '../main.nf'

params.fasta  = null
params.gtf    = null
params.outdir = null

workflow {
    genome_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true), file(params.gtf, checkIfExists: true)])

    empty_ch = Channel.value([[:], []])

    SIMPLEAF_INDEX(genome_ch, empty_ch, empty_ch, empty_ch)
}
