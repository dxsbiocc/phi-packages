#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/modules/esmfold2.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { ESMFOLD2 } from '../main.nf'

params.fasta = null
params.sequences_tsv = null
params.model_dir = null
params.num_recycles = null
params.chunk_size = null
params.args = null
params.outdir = null

workflow {
    def raw_fa = params.fasta ? file(params.fasta, checkIfExists: true) : []
    def fa_list = raw_fa instanceof List ? raw_fa : (raw_fa ? [raw_fa] : [])
    def seq_tsv = params.sequences_tsv ? file(params.sequences_tsv, checkIfExists: true) : []
    def model_dir = params.model_dir ? file(params.model_dir, checkIfExists: true) : []

    ESMFOLD2(
        Channel.value([[id: 'esmfold2'], fa_list, seq_tsv]),
        model_dir,
        params.num_recycles,
        params.chunk_size
    )
}
