#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Composes simpleaf/index directly (vendored module, not a wrapper) so
// callers only supply a genome FASTA/GTF instead of a pre-built index.
nextflow.enable.dsl = 2

include { SIMPLEAF_INDEX } from '../../index/main.nf'
include { SIMPLEAF_QUANT } from '../main.nf'

params.fasta      = null
params.gtf        = null
params.reads_1    = null
params.reads_2    = null
params.chemistry  = '10xv3'
params.resolution = 'cr-like'
params.outdir     = null

workflow {
    genome_ch = Channel.value([[id: 'human'], file(params.fasta, checkIfExists: true), file(params.gtf, checkIfExists: true)])
    empty_ch  = Channel.value([[:], []])

    SIMPLEAF_INDEX(genome_ch, empty_ch, empty_ch, empty_ch)

    reads_ch = Channel.of([
        [id: 'test_10x', single_end: false],
        params.chemistry,
        [file(params.reads_1, checkIfExists: true), file(params.reads_2, checkIfExists: true)]
    ])

    index_ch = SIMPLEAF_INDEX.out.index.combine(SIMPLEAF_INDEX.out.t2g, by: 0)

    cell_filter_ch = Channel.value([[:], 'knee', [], []])
    map_dir_ch     = Channel.value([[:], []])

    SIMPLEAF_QUANT(reads_ch, index_ch, cell_filter_ch, Channel.value(params.resolution), map_dir_ch)
}
