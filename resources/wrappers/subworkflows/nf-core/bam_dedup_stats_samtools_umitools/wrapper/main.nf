#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1. Includes
// the subworkflow only — never another wrapper (section 6).
nextflow.enable.dsl = 2

include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS } from '../main.nf'

params.bam              = null
params.bai              = null
params.single_end       = false
params.get_dedup_stats  = false
params.primary_only     = false
params.outdir           = null

workflow {
    bam_bai_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName, single_end: params.single_end], bam] }
        .combine(Channel.fromPath(params.bai, checkIfExists: true))
        .map { meta, bam, bai -> [meta, bam, bai] }

    BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS(bam_bai_ch, params.get_dedup_stats, params.primary_only)
}
