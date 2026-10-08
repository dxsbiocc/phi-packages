#!/usr/bin/env nextflow
// This file is ../main.nf verbatim, with its 7 relative `include` paths
// adjusted from `./...` to `../...` since this copy lives one directory
// deeper (wrapper/), so the fixed smoke command from
// docs/design/phi-wrapper-agent-composition-design.md section 5 —
// `nextflow run wrapper/main.nf -params-file wrapper/params.json` — works.
// A symlink here instead of a copy does NOT work: Nextflow resolves
// relative `include` paths (and the global `projectDir`) against the
// script path it was literally invoked with, not the symlink's target —
// confirmed by trying it and hitting "Can't find a matching module file
// for include: ./workflows/rnaseq" (the original, unadjusted path).
// Nothing else is changed here: this is the real, unmodified nf-core/rnaseq
// pipeline logic (PIPELINE_INITIALISATION, PREPARE_GENOME, RNASEQ,
// PIPELINE_COMPLETION), minus the Sentieon/Parabricks branches removed
// throughout align_star/quantify_rsem/prepare_genome (see those files).
// ../workflows/rnaseq.nf needed the same kind of `../` fix for its own
// `projectDir`-relative asset paths — see that file's own top-of-file
// comment. There is no `bin/` here (or in ../bin): the one script a
// process actually called by bare name, deseq2_qc.r, now lives inside
// the shared modules/local/differential-expression/qc module itself as a
// Nextflow template instead of `<projectDir>/bin`. The other two
// scripts formerly here (fastq_dir_to_samplesheet.py,
// mqc_features_stat.py) were dead weight — the former a user-facing CLI
// helper no process ever called, the latter fully superseded by
// modules/nf-core/custom/multiqccustombiotype's own templates/ copy —
// and were deleted rather than kept around unused.
// IMPORTANT: this file is a generated derivative of ../main.nf, not a
// live reference — if ../main.nf changes, regenerate this file (sed the
// include paths, re-add this header) rather than hand-editing both copies.
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnaseq
    Website: https://nf-co.re/rnaseq
    Slack  : https://nfcore.slack.com/channels/rnaseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta            = getGenomeAttribute('fasta')
params.additional_fasta = getGenomeAttribute('additional_fasta')
params.transcript_fasta = getGenomeAttribute('transcript_fasta')
params.gff              = getGenomeAttribute('gff')
params.gtf              = getGenomeAttribute('gtf')
params.gene_bed         = getGenomeAttribute('bed12')
params.bbsplit_index    = getGenomeAttribute('bbsplit')
params.sortmerna_index  = getGenomeAttribute('sortmerna')
params.star_index       = getGenomeAttribute('star')
params.rsem_index       = getGenomeAttribute('rsem')
params.hisat2_index     = getGenomeAttribute('hisat2')
params.salmon_index     = getGenomeAttribute('salmon')
params.kallisto_index   = getGenomeAttribute('kallisto')
params.bowtie2_index    = getGenomeAttribute('bowtie2')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQ                  } from '../workflows/rnaseq'
include { PREPARE_GENOME          } from '../subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PIPELINE_COMPLETION     } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { checkMaxContigSize      } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { defineQcTools           } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { isStarIndexLegacy       } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow NFCORE_RNASEQ {

    main:

    //
    // SUBWORKFLOW: Prepare reference genome files
    //
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.gene_bed,
        params.splicesites,
        params.bbsplit_fasta_list,
        params.ribo_database_manifest,
        params.star_index,
        params.rsem_index,
        params.salmon_index,
        params.kallisto_index,
        params.hisat2_index,
        params.bowtie2_index,
        params.bbsplit_index,
        params.sortmerna_index,
        params.kraken_db,
        params.gencode,
        params.gffread_transcript_fasta,
        params.featurecounts_group_type,
        params.aligner,
        params.pseudo_aligner,
        params.skip_gtf_filter,
        params.skip_bbsplit,
        params.remove_ribo_rna ? params.ribo_removal_tool : null,
        params.skip_alignment,
        params.skip_pseudo_alignment,
        params.contaminant_screening,
        params.prokaryotic ?: false,
        isStarIndexLegacy() ?: false
    )

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { fai -> checkMaxContigSize(fai) }
    }

    //
    // WORKFLOW: Run nf-core/rnaseq workflow
    //
    ch_samplesheet = channel.value(file(params.input, checkIfExists: true))

    // Bowtie2 rRNA index is built on-demand inside the fastq_remove_rrna subworkflow
    // rather than in PREPARE_GENOME, to avoid duplicating the rRNA FASTA preparation logic
    ch_bowtie2_rrna_index = channel.empty()

    def qc_tools = defineQcTools(params)

    RNASEQ (
        ch_samplesheet,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.gene_bed,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.rsem_index,
        PREPARE_GENOME.out.hisat2_index,
        PREPARE_GENOME.out.bowtie2_index,
        PREPARE_GENOME.out.salmon_index,
        PREPARE_GENOME.out.kallisto_index,
        PREPARE_GENOME.out.bbsplit_index,
        PREPARE_GENOME.out.rrna_fastas,
        PREPARE_GENOME.out.sortmerna_index,
        ch_bowtie2_rrna_index,
        PREPARE_GENOME.out.splicesites,
        PREPARE_GENOME.out.kraken_db,
        qc_tools
    )

    emit:
    trim_status    = RNASEQ.out.trim_status    // channel: [id, boolean]
    map_status     = RNASEQ.out.map_status     // channel: [id, boolean]
    strand_status  = RNASEQ.out.strand_status  // channel: [id, boolean]
    multiqc_report = RNASEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RNASEQ ()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_RNASEQ.out.multiqc_report,
        NFCORE_RNASEQ.out.trim_status,
        NFCORE_RNASEQ.out.map_status,
        NFCORE_RNASEQ.out.strand_status
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
