#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1. No
// long-read BAM and no explicit mode flags: plain short-read-only
// reference-guided assembly, matching the module's own defaults when
// those are left empty.
nextflow.enable.dsl = 2

include { STRINGTIE_STRINGTIE } from '../main.nf'

params.bam    = null
params.gtf    = null
params.outdir = null

workflow {
    Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.baseName], bam, []] }
        .set { bam_ch }

    STRINGTIE_STRINGTIE(bam_ch, [], file(params.gtf, checkIfExists: true))
}
