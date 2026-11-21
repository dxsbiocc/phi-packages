#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { SAMPLESHEETPARSER_MERGE } from '../main.nf'

params.samplesheets = null
params.target_version = null
params.args = null
params.outdir = null

workflow {
    sheets_ch = Channel.fromPath("${params.samplesheets}/*.csv").collect()

    SAMPLESHEETPARSER_MERGE(sheets_ch.map { s -> [[id: 'merge'], s] }, params.target_version)
}
