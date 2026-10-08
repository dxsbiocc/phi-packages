#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1. No
// existing .fai is ever passed in (samtools faidx builds one from scratch)
// and get_sizes is left off, matching the module's own defaults.
nextflow.enable.dsl = 2

include { SAMTOOLS_FAIDX } from '../main.nf'

params.fasta  = null
params.outdir = null

workflow {
    Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> [[id: fasta.baseName], fasta, []] }
        .set { fasta_ch }

    SAMTOOLS_FAIDX(fasta_ch, false)
}
