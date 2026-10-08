#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// Composes the vendored diamond/makedb module (../../makedb/main.nf) so
// the agent can DIAMOND-BLAST a query against a subject protein FASTA in
// one step, without a separately-built database. See
// docs/design/phi-wrapper-agent-composition-design.md sections 1 and 3.
nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../makedb/main.nf'
include { DIAMOND_BLASTP } from '../main.nf'

params.query   = null
params.subject = null
params.outdir  = null

workflow {
    subject_ch = Channel.value([[id: 'diamonddb'], file(params.subject, checkIfExists: true)])
    empty_ch   = Channel.value([])

    DIAMOND_MAKEDB(subject_ch, empty_ch, empty_ch, empty_ch)

    query_ch = Channel.value([[id: 'query'], file(params.query, checkIfExists: true)])
    db_ch    = DIAMOND_MAKEDB.out.db.map { meta, db -> [[id: 'diamonddb'], db] }

    DIAMOND_BLASTP(query_ch, db_ch, Channel.value(6), Channel.value([]))
}
