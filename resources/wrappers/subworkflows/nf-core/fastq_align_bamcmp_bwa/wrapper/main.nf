#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_ALIGN_BAMCMP_BWA } from '../main.nf'

params.reads = null
params.primary_index = null
params.contaminant_index = null
params.fasta = null
params.fai = null
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    primary_ch = Channel.value([[id: 'primary'], [file(params.primary_index, checkIfExists: true)]])
    contaminant_ch = Channel.value([[id: 'contaminant'], [file(params.contaminant_index, checkIfExists: true)]])
    fasta_fai_ch = params.fasta ? Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)]) : Channel.empty()

    FASTQ_ALIGN_BAMCMP_BWA(reads_ch, primary_ch, contaminant_ch, fasta_fai_ch)
}
