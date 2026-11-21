#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_ALIGN_DEDUP_BWAMEM } from '../main.nf'

params.reads = null
params.fasta = null
params.fai = null
params.bwamem_index = null
params.skip_deduplication = false
params.use_gpu = false
params.output_fmt = "bam"
params.interval_file = null
params.known_sites = null
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    fasta_fai_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])
    index_ch = Channel.value([[id: 'reference'], [file(params.bwamem_index, checkIfExists: true)]])
    interval_ch = params.interval_file ? Channel.value([[id: 'reference'], [file(params.interval_file, checkIfExists: true)]]) : Channel.empty()
    known_ch = params.known_sites ? Channel.value([[id: 'reference'], [file(params.known_sites, checkIfExists: true)]]) : Channel.empty()

    FASTQ_ALIGN_DEDUP_BWAMEM(reads_ch, fasta_fai_ch, index_ch, params.skip_deduplication, params.use_gpu, params.output_fmt, interval_ch, known_ch)
}
