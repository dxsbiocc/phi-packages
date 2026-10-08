#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1. Unlike
// flagstat, `samtools idxstats` needs its .bai staged alongside the BAM
// (by convention, not by name on the command line) — derived by appending
// `.bai` to the raw `params.bam` *string* first, matching how nf-core's own
// test data ships companion index files. Appending it to the already
// `file()`-resolved object instead (`"${bam}.bai"`) silently mangles a
// remote http(s) URL's scheme — confirmed by a real run failing with
// "No such file or directory: /nf-core/test-datasets/...bam.bai" (the
// `https://raw.githubusercontent.com` prefix vanished).
nextflow.enable.dsl = 2

include { SAMTOOLS_IDXSTATS } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    ch_bai = file("${params.bam}.bai", checkIfExists: true)

    Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.baseName], bam, ch_bai] }
        .set { bam_ch }

    SAMTOOLS_IDXSTATS(bam_ch)
}
