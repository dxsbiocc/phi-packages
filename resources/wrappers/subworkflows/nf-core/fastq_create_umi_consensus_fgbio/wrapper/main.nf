#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO } from '../main.nf'

params.reads = null
params.fasta = null
params.fai = null
params.dict = null
params.bwa_index = null
params.groupreadsbyumi_strategy = "Adjacency"
params.aligner = "bwa-mem"
params.duplex = false
params.min_reads = "1"
params.min_baseq = 10
params.max_base_error_rate = 0.025
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    ref_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true), file(params.dict, checkIfExists: true)])
    bwa_ch = params.bwa_index ? Channel.value([[id: 'reference'], file(params.bwa_index, checkIfExists: true)]) : Channel.empty()

    FASTQ_CREATE_UMI_CONSENSUS_FGBIO(reads_ch, ref_ch, bwa_ch, params.groupreadsbyumi_strategy, params.aligner, params.duplex, params.min_reads, params.min_baseq, params.max_base_error_rate)
}
