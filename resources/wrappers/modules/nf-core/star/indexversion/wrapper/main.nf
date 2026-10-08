#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { STAR_INDEXVERSION } from '../main.nf'

params.outdir = null

workflow {
    STAR_INDEXVERSION()
}
