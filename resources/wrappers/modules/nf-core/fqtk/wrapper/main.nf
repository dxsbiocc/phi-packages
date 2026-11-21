#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { FQTK } from '../main.nf'

params.sample_sheet = null
params.fastq_folder = null
params.fastq_names = null
params.read_structures = null
params.args = null
params.outdir = null

workflow {
    def names   = params.fastq_names.tokenize(',')*.trim()
    def structs = params.read_structures.tokenize(',')*.trim()
    if (names.size() != structs.size()) {
        error "fastq_names (${names.size()}) and read_structures (${structs.size()}) must have the same number of comma-separated entries."
    }
    def pairs = [names, structs].transpose()

    FQTK(Channel.value([
        [id: 'fqtk'],
        file(params.sample_sheet, checkIfExists: true),
        file(params.fastq_folder, checkIfExists: true),
        pairs
    ]))
}
