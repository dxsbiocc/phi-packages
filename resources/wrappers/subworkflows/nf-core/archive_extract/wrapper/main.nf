#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { ARCHIVE_EXTRACT } from '../main.nf'

params.archive = null
params.outdir = "results"

workflow {
    archive_ch = Channel
        .fromPath(params.archive, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f] }

    ARCHIVE_EXTRACT(archive_ch)
}
