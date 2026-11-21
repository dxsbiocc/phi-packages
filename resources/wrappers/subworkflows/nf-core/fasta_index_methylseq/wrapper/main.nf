#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTA_INDEX_METHYLSEQ } from '../main.nf'

params.fasta = null
params.fasta_index = null
params.bismark_index = null
params.bwameth_index = null
params.bwamem_index = null
params.aligner = "bismark"
params.collecthsmetrics = false
params.use_mem2 = false
params.outdir = "results"

workflow {
    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])
    fai_ch = params.fasta_index ? Channel.value([[id: 'reference'], file(params.fasta_index, checkIfExists: true)]) : []
    bismark_ch = params.bismark_index ? Channel.value([[id: 'reference'], file(params.bismark_index, checkIfExists: true)]) : []
    bwameth_ch = params.bwameth_index ? Channel.value([[id: 'reference'], file(params.bwameth_index, checkIfExists: true)]) : []
    bwamem_ch = params.bwamem_index ? Channel.value([[id: 'reference'], file(params.bwamem_index, checkIfExists: true)]) : []

    FASTA_INDEX_METHYLSEQ(fasta_ch, fai_ch, bismark_ch, bwameth_ch, bwamem_ch, params.aligner, params.collecthsmetrics, params.use_mem2)
}
