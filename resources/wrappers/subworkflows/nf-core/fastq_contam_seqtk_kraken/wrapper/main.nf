#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_CONTAM_SEQTK_KRAKEN } from '../main.nf'

params.reads = null
params.sample_size = "10000"
params.kraken2_db = null
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    kraken2_db_ch = Channel.value(file(params.kraken2_db, checkIfExists: true))

    FASTQ_CONTAM_SEQTK_KRAKEN(reads_ch, params.sample_size, kraken2_db_ch)
}
