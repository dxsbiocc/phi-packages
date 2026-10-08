#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// Composes the vendored blast/makeblastdb module (../../makeblastdb/main.nf)
// so the agent can BLAST a protein query against a raw subject protein
// FASTA in one step, without a separately-built database. See
// docs/design/phi-wrapper-agent-composition-design.md sections 1 and 3.
nextflow.enable.dsl = 2

include { BLAST_MAKEBLASTDB } from '../../makeblastdb/main.nf'
include { BLAST_BLASTP }      from '../main.nf'

params.query   = null
params.subject = null
params.outdir  = null

workflow {
    subject_ch   = Channel.value([[id: 'blastdb'], file(params.subject, checkIfExists: true)])
    taxid_map_ch = Channel.value([])

    BLAST_MAKEBLASTDB(subject_ch, taxid_map_ch)

    query_ch = Channel.value([[id: 'query'], file(params.query, checkIfExists: true)])
    db_ch    = BLAST_MAKEBLASTDB.out.db.map { meta, db -> [[id: 'blastdb'], db] }

    BLAST_BLASTP(query_ch, db_ch, 'tsv')
}
