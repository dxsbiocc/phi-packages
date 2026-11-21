#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTA_VCLUST_PREFILTER_ALIGN_CLUSTER } from '../main.nf'

params.fasta = null
params.save_alignment = false
params.metric = "ani"
params.tani = 0.95
params.gani = 0.95
params.ani = 0.95
params.outdir = "results"

workflow {
    fasta_ch = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .collect()
        .map { files -> [[id: 'vclust'], files] }

    FASTA_VCLUST_PREFILTER_ALIGN_CLUSTER(fasta_ch, params.save_alignment, params.metric, params.tani, params.gani, params.ani)
}
