nextflow.enable.dsl = 2

include { FOLDSEEK_EASYSEARCH } from '../main.nf'

params.structures = null
params.target_db = null
params.evalue = null
params.novelty_tm_cutoff = null
params.args = null
params.outdir = null

workflow {
    if (!params.structures) {
        error "Missing required parameter: --structures"
    }
    if (!params.target_db) {
        error "Missing required parameter: --target_db"
    }

    def struct_files = files(params.structures)
    if (!struct_files) {
        error "No structure files matched --structures: ${params.structures}"
    }

    ch_structures = Channel.value([[id: 'designs'], struct_files])
    ch_target_db = Channel.value(file(params.target_db, checkIfExists: true))

    FOLDSEEK_EASYSEARCH(
        ch_structures,
        ch_target_db,
        params.evalue,
        params.novelty_tm_cutoff
    )
}
