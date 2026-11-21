#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { FASTQSCREEN_BUILDFROMINDEX } from '../main.nf'

params.indexes = null
params.outdir = null

workflow {
    dirs_ch = Channel.fromPath(params.indexes, checkIfExists: true)
        .map { d -> d.listFiles().findAll { it.isDirectory() }.toList() }

    FASTQSCREEN_BUILDFROMINDEX(dirs_ch.map { ds -> ds*.name }, dirs_ch)
}
