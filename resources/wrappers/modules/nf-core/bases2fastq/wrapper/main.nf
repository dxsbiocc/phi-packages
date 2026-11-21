#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { BASES2FASTQ } from '../main.nf'

params.run_manifest = null
params.run_dir = null
params.args = null
params.outdir = null

workflow {
    BASES2FASTQ(Channel.value([
        [id: 'run'],
        file(params.run_manifest, checkIfExists: true),
        file(params.run_dir, checkIfExists: true)
    ]))
}
