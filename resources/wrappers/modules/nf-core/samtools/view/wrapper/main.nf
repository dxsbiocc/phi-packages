#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1. No
// reference/qname/bed filters and no re-indexing: passes the BAM through
// unfiltered, matching the module's own defaults when those are left
// empty. The id gets a "_view" suffix — the module errors if its output
// filename would collide with the input's, which a bare `baseName` id
// does for a same-format passthrough.
nextflow.enable.dsl = 2

include { SAMTOOLS_VIEW } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: "${bam.baseName}_view"], bam, []] }
        .set { bam_ch }

    ch_fasta = channel.value([[:], [], []])
    ch_qname = channel.value([[:], []])
    ch_bed   = channel.value([[:], []])

    SAMTOOLS_VIEW(bam_ch, ch_fasta, ch_qname, ch_bed, '')
}
