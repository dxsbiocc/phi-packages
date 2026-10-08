#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1. No
// genome FASTA (decoy sequences): a transcript-only index, matching the
// module's own default when that input is left empty.
nextflow.enable.dsl = 2

include { SALMON_INDEX } from '../main.nf'

params.transcript_fasta = null
params.outdir           = null

workflow {
    SALMON_INDEX(
        [],
        file(params.transcript_fasta, checkIfExists: true)
    )
}
