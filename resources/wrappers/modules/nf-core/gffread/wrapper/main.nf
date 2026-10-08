#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GFFREAD } from '../main.nf'

params.gff    = null
params.outdir = null

workflow {
    gff_ch = Channel
        .fromPath(params.gff, checkIfExists: true)
        .map { gff -> [[id: gff.simpleName], gff] }

    GFFREAD(gff_ch, [])
}
