nextflow.enable.dsl = 2

include { FOLDSEEK_EASYCLUSTER } from '../main.nf'

params.structures = null
params.coverage = null
params.tmscore_threshold = null
params.args = null
params.outdir = null

workflow {
    if (!params.structures) {
        error "Missing required parameter: --structures"
    }

    def struct_files = files(params.structures)
    if (!struct_files) {
        error "No structure files matched --structures: ${params.structures}"
    }

    ch_structures = Channel.value([[id: 'designs'], struct_files])

    FOLDSEEK_EASYCLUSTER(
        ch_structures,
        params.coverage,
        params.tmscore_threshold
    )
}
