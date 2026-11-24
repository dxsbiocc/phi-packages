#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/modules/foldseek.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { FOLDSEEK } from '../main.nf'

params.structures = null
params.target_db = null
params.mode = null
params.evalue = null
params.coverage = null
params.tmscore_threshold = null
params.novelty_tm_cutoff = null
params.args = null
params.outdir = null

workflow {
    if (params.mode == 'search' && !params.target_db) {
        error "foldseek requires `target_db` when `mode` is set to 'search'."
    }
    def struct_files = file(params.structures, checkIfExists: true)
    def struct_list = struct_files instanceof List ? struct_files : [struct_files]
    def db = params.target_db ? file(params.target_db, checkIfExists: true) : []

    FOLDSEEK(
        Channel.value([[id: 'foldseek'], struct_list]),
        db,
        params.mode ?: 'cluster',
        params.evalue,
        params.coverage,
        params.tmscore_threshold,
        params.novelty_tm_cutoff
    )
}
