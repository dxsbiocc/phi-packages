#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/kneaddata, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { KNEADDATA } from '../main.nf'

params.reads = null
params.db = null
params.args = null
params.outdir = null

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, size: -1)
        .map { sample, reads -> [[id: sample, single_end: reads instanceof List ? reads.size() == 1 : true], reads] }

    KNEADDATA(reads_ch, file(params.db, checkIfExists: true))
}
