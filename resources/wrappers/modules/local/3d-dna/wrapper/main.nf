#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/3d-dna, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { THREEDDNA } from '../main.nf'

params.fasta = null
params.mnd = null
params.mode = null
params.input_size = null
params.rounds = null
params.stage = null
params.args = null
params.outdir = null

workflow {
    THREEDDNA(Channel.value([
        [id: 'threeddna'],
        file(params.fasta, checkIfExists: true),
        file(params.mnd, checkIfExists: true)
    ]))
}
