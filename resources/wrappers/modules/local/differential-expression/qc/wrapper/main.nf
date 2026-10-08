#!/usr/bin/env nextflow
// Thin agent-facing adapter over the DESeq2 QC module at ../main.nf —
// sample relationships from a counts matrix: PCA and hierarchical-
// clustering plots, no contrast/hypothesis testing. See
// docs/design/phi-wrapper-agent-composition-design.md section 1.
//
// For testing a specific contrast, use ../../deseq2/wrapper/ (or limma/,
// timeseries/) in this same differential-expression family.
nextflow.enable.dsl = 2

include { DESEQ2_QC } from '../main.nf'

params.counts = null
params.outdir = null

workflow {
    counts_ch = Channel.fromPath(params.counts, checkIfExists: true)
    // The shared process also serves nf-core/rnaseq, which supplies real
    // MultiQC headers. This direct adapter does not
    // publish those optional fragments, so one existing small file supplies
    // the two required path inputs.
    support_file_ch = Channel.value(file("${projectDir}/params.json", checkIfExists: true))

    DESEQ2_QC(counts_ch, support_file_ch, support_file_ch)
}
