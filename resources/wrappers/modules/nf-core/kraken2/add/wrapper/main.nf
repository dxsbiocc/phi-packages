#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { KRAKEN2_ADD } from '../main.nf'
include { GUNZIP } from '../../../gunzip/main.nf'

params.fasta      = null
params.proteome   = null
params.names      = null
params.nodes      = null
params.accession2taxid = null
params.outdir     = null

workflow {
    GUNZIP(Channel.value([[], file(params.accession2taxid, checkIfExists: true)]))
    KRAKEN2_ADD(Channel.value([[id: 'test'], [file(params.fasta, checkIfExists: true), file(params.proteome, checkIfExists: true)]]),
                file(params.names, checkIfExists: true), file(params.nodes, checkIfExists: true),
                GUNZIP.out.gunzip.map { it[1] }, [])
}
