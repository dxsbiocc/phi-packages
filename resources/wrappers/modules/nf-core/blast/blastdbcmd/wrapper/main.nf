#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BLAST_BLASTDBCMD } from '../main.nf'
include { BLAST_MAKEBLASTDB } from '../../makeblastdb/main.nf'

params.fasta      = null
params.entry      = 'ENSSASP00005000002.1'
params.outdir     = null

workflow {
    BLAST_MAKEBLASTDB(Channel.value([[id: 'blastdb'], file(params.fasta, checkIfExists: true)]), Channel.value([]))
    BLAST_BLASTDBCMD(Channel.value([[id: 'result'], params.entry, []]), BLAST_MAKEBLASTDB.out.db)
}
