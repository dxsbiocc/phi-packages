#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS } from '../main.nf'

params.reads = null
params.db = null
params.taxid = "9606"
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    db_ch = Channel.value(file(params.db, checkIfExists: true))

    FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS(reads_ch, db_ch, params.taxid)
}
