#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_ALIGN_CHROMAP } from '../main.nf'

params.reads = null
params.index = null
params.fasta = null
params.fai = null
params.barcodes = null
params.whitelist = null
params.chr_order = null
params.pairs_chr_order = null
params.update_readgroups = false
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    index_ch = Channel.value([[id: 'reference'], [file(params.index, checkIfExists: true)]])
    fasta_fai_ch = Channel.value([[id: 'reference'], [file(params.fasta, checkIfExists: true)], [file(params.fai, checkIfExists: true)]])
    barcodes_ch = params.barcodes ? Channel.value(file(params.barcodes, checkIfExists: true)) : Channel.empty()
    whitelist_ch = params.whitelist ? Channel.value(file(params.whitelist, checkIfExists: true)) : Channel.empty()
    chr_order_ch = params.chr_order ? Channel.value(file(params.chr_order, checkIfExists: true)) : Channel.empty()
    pairs_ch = params.pairs_chr_order ? Channel.value(file(params.pairs_chr_order, checkIfExists: true)) : Channel.empty()

    FASTQ_ALIGN_CHROMAP(reads_ch, index_ch, fasta_fai_ch, barcodes_ch, whitelist_ch, chr_order_ch, pairs_ch, params.update_readgroups)
}
