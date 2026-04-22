//
// FORMATTING SUB-WORKFLOW
//

include { FORMATTER } from "../../../subworkflows/local/formatter/main"
include { PREPARE_CLUSTER as PREPARE_CLUSTER_PYCLONE } from '../../../modules/local/prepare_cluster/main'
include { PREPARE_CLUSTER as PREPARE_CLUSTER_MOBSTER } from '../../../modules/local/prepare_cluster/main'
include { PREPARE_CLUSTER as PREPARE_CLUSTER_VIBER } from '../../../modules/local/prepare_cluster/main'
include { ASSIGN_CLUSTER as ASSIGN_CLUSTER_VIBER } from '../../../modules/local/assign_cluster/main'
include { ASSIGN_CLUSTER as ASSIGN_CLUSTER_MOBSTER } from '../../../modules/local/assign_cluster/main'
include { ASSIGN_CLUSTER as ASSIGN_CLUSTER_PYCLONE } from '../../../modules/local/assign_cluster/main'

workflow ASSIGN_SIGNATURE {
    take:
        pyclone_fit
        rds_join
        mobster_fit
        viber_fit
        sigprofiler_fit

    main:
        ch_versions = Channel.empty()
        table_pyclone = null
        table_mobster = null
        table_viber = null
        assign_pyclone = null
        assign_mobster = null
        assign_viber = null

        if (params.download_sigprofiler_genome) {
            genome_path = channel.fromPath('opt/null')
        } else {
            genome_path = params.genome_installed_path
        }

        if (params.tools && params.tools.split(",").contains("mobster")) {
            mobster_combined = mobster_fit.map { meta, fit -> tuple(meta, fit, []) }
            PREPARE_CLUSTER_MOBSTER(mobster_combined)
            table_mobster = PREPARE_CLUSTER_MOBSTER.out.signature_table
            ch_versions = ch_versions.mix(PREPARE_CLUSTER_MOBSTER.out.versions)

            input = table_mobster.combine(sigprofiler_fit).combine(genome_path)
            ASSIGN_CLUSTER_MOBSTER(input, 'mobster')
            assign_mobster = ASSIGN_CLUSTER_MOBSTER.out.results_sigprofiler
        }

        if (params.tools && params.tools.split(",").contains("viber")) {
            viber_combined = viber_fit.map { meta, fit -> tuple(meta, fit, []) }
            PREPARE_CLUSTER_VIBER(viber_combined)
            table_viber = PREPARE_CLUSTER_VIBER.out.signature_table
            ch_versions = ch_versions.mix(PREPARE_CLUSTER_VIBER.out.versions)

            input = table_viber.combine(sigprofiler_fit).combine(genome_path)
            ASSIGN_CLUSTER_VIBER(input, 'viber')
        }

        if (params.tools && params.tools.split(",").contains("pyclone-vi")) {
            FORMATTER(rds_join, "rds")
            ch_versions = ch_versions.mix(FORMATTER.out.versions)

            pyclone_combined = pyclone_fit.join(FORMATTER.out.out_data, by: 0).map {meta, fit, data, samples -> tuple(meta, fit, data)}
            PREPARE_CLUSTER_PYCLONE(pyclone_combined)
            table_pyclone = PREPARE_CLUSTER_PYCLONE.out.signature_table
            ch_versions = ch_versions.mix(PREPARE_CLUSTER_PYCLONE.out.versions)

            input = table_pyclone.combine(sigprofiler_fit).combine(genome_path)
            ASSIGN_CLUSTER_PYCLONE(input, 'pyclonevi')
            assign_pyclone = ASSIGN_CLUSTER_PYCLONE.out.results_sigprofiler
        }

    emit:
        table_pyclone
        table_viber
        table_mobster
        assign_pyclone
        assign_viber
        assign_mobster
        ch_versions
}
