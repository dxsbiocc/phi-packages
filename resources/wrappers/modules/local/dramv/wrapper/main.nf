#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/dramv, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { DRAMV } from '../main.nf'

params.fasta = null
params.affi = null
params.config = null
params.distill = null
params.args = null
params.outdir = null

// Requires a DRAM installation with databases configured (DRAM-setup.py
// prepare_databases) and a config file whose amg_database is a real path.
workflow {
    DRAMV(
        Channel.value([[id: 'dramv'],
                       file(params.fasta, checkIfExists: true),
                       file(params.affi, checkIfExists: true)]),
        params.config ? file(params.config, checkIfExists: true) : [],
        params.distill
    )
}
