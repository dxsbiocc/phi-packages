#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1. Includes
// the subworkflow only — never another wrapper (section 6).
nextflow.enable.dsl = 2

include { BAM_STRINGTIE_MERGE } from '../main.nf'

params.bams         = null
params.gtf          = null
params.strandedness = 'reverse'
params.outdir       = null

workflow {
    bam_ch = Channel
        .fromPath(params.bams, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName, strandedness: params.strandedness], bam] }

    gtf_ch = Channel
        .value([[id: 'annotation'], file(params.gtf, checkIfExists: true)])

    BAM_STRINGTIE_MERGE(bam_ch, gtf_ch)
}
