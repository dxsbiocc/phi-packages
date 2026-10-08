#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// Composes the vendored mmseqs/createdb module (../../createdb/main.nf)
// twice so the agent can search a raw query FASTA against a raw target
// FASTA in one step, without separately-built databases. See
// docs/design/phi-wrapper-agent-composition-design.md sections 1 and 3.
nextflow.enable.dsl = 2

include { MMSEQS_CREATEDB as MMSEQS_CREATEDB_QUERY  } from '../../createdb/main.nf'
include { MMSEQS_CREATEDB as MMSEQS_CREATEDB_TARGET } from '../../createdb/main.nf'
include { MMSEQS_SEARCH }                              from '../main.nf'

params.query   = null
params.target  = null
params.outdir  = null

workflow {
    query_ch  = Channel.value([[id: 'querydb', single_end: false], file(params.query, checkIfExists: true)])
    target_ch = Channel.value([[id: 'targetdb', single_end: false], file(params.target, checkIfExists: true)])

    MMSEQS_CREATEDB_QUERY(query_ch)
    MMSEQS_CREATEDB_TARGET(target_ch)

    MMSEQS_SEARCH(MMSEQS_CREATEDB_QUERY.out.db, MMSEQS_CREATEDB_TARGET.out.db)
}
