#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_ADDMOSTSEVERECONSEQUENCE } from '../main.nf'


params.vcf        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], file(params.vcf, checkIfExists: true)])
    csq_ch = Channel.of(['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification', 'feature_elongation', 'feature_truncation', 'inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant', 'splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant', 'splice_polypyrimidine_tract_variant', 'incomplete_terminal_codon_variant', 'start_retained_variant', 'stop_retained_variant', 'synonymous_variant', 'coding_sequence_variant', 'mature_miRNA_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', 'non_coding_transcript_exon_variant', 'intron_variant', 'NMD_transcript_variant', 'non_coding_transcript_variant', 'coding_transcript_variant', 'upstream_gene_variant', 'downstream_gene_variant', 'TFBS_ablation', 'TFBS_amplification', 'TF_binding_site_variant', 'regulatory_region_ablation', 'regulatory_region_amplification', 'regulatory_region_variant', 'intergenic_variant', 'sequence_variant'])
        .collectFile(name: 'variant_consequences.txt', newLine: true)
        .map { file -> [[id: file.simpleName], file] }

    CUSTOM_ADDMOSTSEVERECONSEQUENCE(in_ch, csq_ch)
}
