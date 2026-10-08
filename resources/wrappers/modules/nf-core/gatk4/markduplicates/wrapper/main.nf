#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// The `.markdup` id suffix keeps the output BAM name different from the
// input's — same reasoning as other same-name-collision fixes in this tree.
nextflow.enable.dsl = 2

include { GATK4_MARKDUPLICATES } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: "${bam.simpleName}.markdup", single_end: false], bam] }

    GATK4_MARKDUPLICATES(bam_ch, [], [])
}
