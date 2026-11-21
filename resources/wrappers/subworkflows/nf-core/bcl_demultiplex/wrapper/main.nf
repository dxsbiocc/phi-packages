#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { BCL_DEMULTIPLEX } from '../main.nf'

params.samplesheet = null
params.flowcell = null
params.run_id = "flowcell"
params.lane = 1
params.demultiplexer = "bclconvert"
params.outdir = "results"

workflow {
    flowcell_ch = Channel.value([[id: params.run_id, lane: params.lane], file(params.samplesheet, checkIfExists: true), file(params.flowcell, checkIfExists: true)])

    BCL_DEMULTIPLEX(flowcell_ch, params.demultiplexer)
}
