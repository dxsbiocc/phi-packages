#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { SAMPLESHEETPARSER_DIFF } from '../main.nf'

params.old_sheet = null
params.new_sheet = null
params.args = null
params.outdir = null

workflow {
    SAMPLESHEETPARSER_DIFF(Channel.value([
        [id: 'diff'],
        file(params.old_sheet, checkIfExists: true),
        file(params.new_sheet, checkIfExists: true)
    ]))
}
