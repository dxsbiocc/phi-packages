#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTA_GXF_BUSCO_PLOT } from '../main.nf'

params.fasta = null
params.gxf = null
params.mode = "genome"
params.lineages = "auto"
params.busco_lineages_path = null
params.busco_config = null
params.busco_cleanup = true
params.outdir = "results"

workflow {
    fasta_ch = Channel.fromPath(params.fasta, checkIfExists: true).map { f -> [[id: 'busco'], f] }
    gxf_ch = params.gxf ? Channel.value([[id: 'busco'], file(params.gxf, checkIfExists: true)]) : Channel.empty()

    FASTA_GXF_BUSCO_PLOT(fasta_ch, gxf_ch, params.mode, params.lineages.tokenize(','), params.busco_lineages_path ?: [], params.busco_config ?: [], params.busco_cleanup)
}
