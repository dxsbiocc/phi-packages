#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTA_INDEX_BISMARK_BWAMETH } from '../main.nf'

params.fasta = null
params.fai = null
params.bismark_index = null
params.bwameth_index = null
params.aligner = "bismark"
params.use_mem2 = false
params.outdir = "results"

workflow {
    fasta_fai_ch = Channel.value([[id: 'reference'], [file(params.fasta, checkIfExists: true)], params.fai ? [file(params.fai, checkIfExists: true)] : []])
    bismark_ch = params.bismark_index ? Channel.value([[id: 'reference'], [file(params.bismark_index, checkIfExists: true)]]) : Channel.empty()
    bwameth_ch = params.bwameth_index ? Channel.value([[id: 'reference'], [file(params.bwameth_index, checkIfExists: true)]]) : Channel.empty()

    FASTA_INDEX_BISMARK_BWAMETH(fasta_fai_ch, bismark_ch, bwameth_ch, params.aligner, params.use_mem2)
}
