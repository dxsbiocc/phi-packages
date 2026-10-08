#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// Composes the vendored mmseqs/createdb module (../../createdb/main.nf)
// so the agent can cluster a raw FASTA in one step, without a
// separately-built database. See
// docs/design/phi-wrapper-agent-composition-design.md sections 1 and 3.
nextflow.enable.dsl = 2

include { MMSEQS_CREATEDB } from '../../createdb/main.nf'
include { MMSEQS_CLUSTER }  from '../main.nf'

params.fasta  = null
params.outdir = null

workflow {
    fasta_ch = Channel.value([[id: 'mmseqsdb', single_end: false], file(params.fasta, checkIfExists: true)])

    MMSEQS_CREATEDB(fasta_ch)

    db_ch = MMSEQS_CREATEDB.out.db.map { meta, db -> [[id: 'mmseqsdb_clustered'], db] }

    MMSEQS_CLUSTER(db_ch)
}
