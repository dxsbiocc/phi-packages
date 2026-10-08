#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1. The bai
// input is never actually read by the module's own script (plain
// `samtools flagstat $bam`), so it's always passed empty here.
nextflow.enable.dsl = 2

include { SAMTOOLS_FLAGSTAT } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.baseName], bam, []] }
        .set { bam_ch }

    SAMTOOLS_FLAGSTAT(bam_ch)
}
