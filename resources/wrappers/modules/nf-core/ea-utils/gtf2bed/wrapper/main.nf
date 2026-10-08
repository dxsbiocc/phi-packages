#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { EAUTILS_GTF2BED } from '../main.nf'

params.gtf    = null
params.outdir = null

workflow {
    gtf_ch = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .map { gtf -> [[id: gtf.simpleName], gtf] }

    EAUTILS_GTF2BED(gtf_ch)
}
