#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// Composes the vendored picard/sortsam module (../../sortsam/main.nf) —
// FilterSamReads asserts queryname order internally, and the module
// itself doesn't sort, so the upstream test queryname-sorts first. See
// docs/design/phi-wrapper-agent-composition-design.md sections 1 and 3.
nextflow.enable.dsl = 2

include { PICARD_SORTSAM }       from '../../sortsam/main.nf'
include { PICARD_FILTERSAMREADS } from '../main.nf'

params.bam    = null
params.filter = 'includeAligned'
params.outdir = null

workflow {
    bam_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true)])

    PICARD_SORTSAM(bam_ch, 'queryname')

    sorted_ch = PICARD_SORTSAM.out.bam.map { meta, bam -> [[id: "${meta.id}_filtered"], bam, []] }
    fasta_ch  = Channel.value([[:], []])

    PICARD_FILTERSAMREADS(sorted_ch, fasta_ch, params.filter)
}
