//
// FORMATTING SUB-WORKFLOW
//

include { FORMATTER } from "../../../subworkflows/local/formatter/main"
include { PREPARE_CLUSTER as PREPARE_CLUSTER_PYCLONE } from '../../../modules/local/prepare_cluster/main'
include { PREPARE_CLUSTER as PREPARE_CLUSTER_MOBSTER } from '../../../modules/local/prepare_cluster/main'
include { PREPARE_CLUSTER as PREPARE_CLUSTER_VIBER } from '../../../modules/local/prepare_cluster/main'


workflow ASSIGN_SIGNATURE {
    take:
        pyclone_fit
        rds_join
        mobster_fit
        viber_fit
        sigprofiler_fit

    main:

        FORMATTER(rds_join, "rds")
        pyclone_combined = pyclone_fit.join(FORMATTER.out.out_data, by: 0).map {meta, fit, data, samples -> tuple(meta, fit, data)}
        viber_combined = viber_fit.map { meta, fit -> tuple(meta, fit, []) }
        mobster_combined = mobster_fit.map { meta, fit -> tuple(meta, fit, []) }

        PREPARE_CLUSTER_PYCLONE(pyclone_combined)
        //PREPARE_CLUSTER_VIBER(viber_combined)
        //PREPARE_CLUSTER_MOBSTER(mobster_combined)

    emit:
        null
}
