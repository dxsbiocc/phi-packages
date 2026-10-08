#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1.
//
// Composes RSEM_PREPAREREFERENCE + RSEM_CALCULATEEXPRESSION internally —
// same reasoning as the other aligner/quantifier wrappers in this tree: no
// pre-built index in test-datasets to just point this wrapper's default
// params.json at.
//
// Two explicit read-file params instead of one `reads` glob: Nextflow
// rejects glob patterns entirely for `https://` sources.
nextflow.enable.dsl = 2

include { RSEM_PREPAREREFERENCE     } from '../../preparereference/main.nf'
include { RSEM_CALCULATEEXPRESSION  } from '../main.nf'

params.reads_1 = null
params.reads_2 = null
params.fasta   = null
params.gtf     = null
params.outdir  = null

workflow {
    RSEM_PREPAREREFERENCE(
        file(params.fasta, checkIfExists: true),
        file(params.gtf, checkIfExists: true)
    )

    reads_ch = channel.of([
        [id: 'test', single_end: false],
        [file(params.reads_1, checkIfExists: true), file(params.reads_2, checkIfExists: true)]
    ])

    RSEM_CALCULATEEXPRESSION(reads_ch, RSEM_PREPAREREFERENCE.out.index)
}
