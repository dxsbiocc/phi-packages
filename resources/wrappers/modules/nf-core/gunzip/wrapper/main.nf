#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GUNZIP } from '../main.nf'

params.archive = null
params.outdir  = null

workflow {
    archive_ch = Channel
        .fromPath(params.archive, checkIfExists: true)
        .map { archive -> [[id: archive.simpleName], archive] }

    GUNZIP(archive_ch)
}
