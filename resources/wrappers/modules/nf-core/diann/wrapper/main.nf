#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { DIANN } from '../main.nf'

params.ms_files = null
params.ms_file_names = null
params.fasta = null
params.library = null
params.quant = null
params.args = null
params.outdir = null

// The vendored module refuses -profile conda; use docker/singularity.
workflow {
    def ms = params.ms_files
        ? file(params.ms_files, checkIfExists: true).listFiles().findAll { it.isFile() }.toList()
        : []
    def names = params.ms_file_names ? params.ms_file_names.tokenize(',')*.trim() : []
    def quant = params.quant
        ? file(params.quant, checkIfExists: true).listFiles().findAll { it.isFile() }.toList()
        : []

    DIANN(Channel.value([
        [id: 'diann'], ms, names,
        params.fasta   ? file(params.fasta, checkIfExists: true)   : [],
        params.library ? file(params.library, checkIfExists: true) : [],
        quant
    ]))
}
