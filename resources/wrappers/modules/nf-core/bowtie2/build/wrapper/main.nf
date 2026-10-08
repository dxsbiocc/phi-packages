#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../main.nf'

params.fasta  = null
params.outdir = null

workflow {
    Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> [[id: fasta.baseName], fasta] }
        .set { fasta_ch }

    BOWTIE2_BUILD(fasta_ch)
}
