#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_DECONTAMINATE_DEACON_HOSTILE } from '../main.nf'

params.reads = null
params.fasta = null
params.reference = null
params.reference_name = "contaminant"
params.index_name = null
params.decontaminator = "hostile"
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    fasta_ch = params.fasta ? Channel.value([[id: params.reference_name], [file(params.fasta, checkIfExists: true)]]) : Channel.empty()
    reference_ch = params.reference ? Channel.value([params.reference_name, file(params.reference, checkIfExists: true)]) : Channel.empty()

    FASTQ_DECONTAMINATE_DEACON_HOSTILE(reads_ch, fasta_ch, reference_ch, params.index_name ?: null, params.decontaminator)
}
