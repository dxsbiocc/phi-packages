#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { BAM_VARIANT_DEMIX_BOOT_FREYJA } from '../main.nf'

params.bam = null
params.fasta = null
params.skip_boot = true
params.repeats = 100
params.db_name = "freyja"
params.barcodes = null
params.lineages_meta = null
params.lineages_topology = null
params.outdir = "results"

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f] }
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    barcodes_ch = params.barcodes ? Channel.value(file(params.barcodes, checkIfExists: true)) : Channel.empty()
    meta_ch = params.lineages_meta ? Channel.value(file(params.lineages_meta, checkIfExists: true)) : Channel.empty()
    topo_ch = params.lineages_topology ? Channel.value(file(params.lineages_topology, checkIfExists: true)) : Channel.empty()

    BAM_VARIANT_DEMIX_BOOT_FREYJA(bam_ch, fasta_ch, params.skip_boot, params.repeats, params.db_name, barcodes_ch, meta_ch, topo_ch)
}
