#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// Composes the vendored mmseqs/createdb module (../../createdb/main.nf)
// so the agent can search a query FASTA against a raw target FASTA in
// one step, without a separately-built database. See
// docs/design/phi-wrapper-agent-composition-design.md sections 1 and 3.
nextflow.enable.dsl = 2

include { MMSEQS_CREATEDB }   from '../../createdb/main.nf'
include { MMSEQS_EASYSEARCH } from '../main.nf'

params.query   = null
params.target  = null
params.outdir  = null

workflow {
    target_ch = Channel.value([[id: 'mmseqsdb', single_end: false], file(params.target, checkIfExists: true)])

    MMSEQS_CREATEDB(target_ch)

    query_ch = Channel.value([[id: 'query', single_end: true], file(params.query, checkIfExists: true)])

    MMSEQS_EASYSEARCH(query_ch, MMSEQS_CREATEDB.out.db)
}
