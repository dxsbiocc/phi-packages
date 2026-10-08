#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BBMAP_FILTERBYNAME } from '../main.nf'

params.reads1 = null
params.reads2 = null
params.names  = ''
params.outdir = null

workflow {
    reads_ch = Channel.value([[id: 'filtered', single_end: false], [file(params.reads1, checkIfExists: true), file(params.reads2, checkIfExists: true)]])

    BBMAP_FILTERBYNAME(reads_ch, params.names, 'fastq.gz', false)
}
