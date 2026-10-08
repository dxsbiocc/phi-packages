#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CAT_FASTQ } from '../main.nf'

params.reads  = null
params.outdir = null

workflow {
    Channel
        .fromPath(params.reads, checkIfExists: true)
        .toSortedList()
        .map { files -> [[id: 'merged', single_end: true], files] }
        .set { reads_ch }

    CAT_FASTQ(reads_ch)
}
