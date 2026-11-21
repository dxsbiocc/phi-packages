#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTA_NEWICK_EPANG_GAPPA } from '../main.nf'

params.queryseqfile = null
params.refseqfile = null
params.refphylogeny = null
params.hmmfile = null
params.alignmethod = "hmmer"
params.model = null
params.outdir = "results"

workflow {
    pp_ch = Channel.value([
        meta: [id: 'epang'],
        data: [
            alignmethod: params.alignmethod,
            queryseqfile: file(params.queryseqfile, checkIfExists: true),
            refseqfile: file(params.refseqfile, checkIfExists: true),
            refphylogeny: file(params.refphylogeny, checkIfExists: true),
            hmmfile: params.hmmfile ? file(params.hmmfile, checkIfExists: true) : [],
            model: params.model ?: null
        ]
    ])

    FASTA_NEWICK_EPANG_GAPPA(pp_ch)
}
