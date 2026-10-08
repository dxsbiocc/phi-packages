#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
//
// Only the "build an index" mode is exposed (only_build_index=true, no
// reads) — see the module's own test suite, which only exercises indexing
// with this tiny data (splitting real reads needs a much larger reference
// set to say anything meaningful). Splitting itself reuses the same
// BBMAP_BBSPLIT process with an existing index, a separate use case.
nextflow.enable.dsl = 2

include { BBMAP_BBSPLIT } from '../main.nf'

params.primary_ref = null
params.other_refs  = null
params.other_names = 'human'
params.outdir      = null

workflow {
    reads_ch = Channel.value([[:], []])
    index_ch = Channel.value([])

    primary_ref_ch = Channel.fromPath(params.primary_ref, checkIfExists: true)

    other_names = params.other_names.split(',') as List
    other_ref_ch = Channel
        .fromPath(params.other_refs, checkIfExists: true)
        .collect()
        .map { refs -> [other_names, refs] }

    BBMAP_BBSPLIT(reads_ch, index_ch, primary_ref_ch, other_ref_ch, true)
}
