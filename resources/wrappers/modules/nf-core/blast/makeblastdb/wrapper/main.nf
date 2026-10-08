#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BLAST_MAKEBLASTDB } from '../main.nf'

params.fasta   = null
params.dbtype  = 'nucl'
params.outdir  = null

workflow {
    fasta_ch     = Channel.value([[id: 'blastdb'], file(params.fasta, checkIfExists: true)])
    taxid_map_ch = Channel.value([])

    BLAST_MAKEBLASTDB(fasta_ch, taxid_map_ch)
}
