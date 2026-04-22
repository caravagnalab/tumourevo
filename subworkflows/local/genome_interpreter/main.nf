//
// GENOME INTERPRETER WORKFLOW
//

include { COHORT_QC } from "../../../modules/local/cohort_qc/main"
include { COHORT_MUTATIONS } from "../../../modules/local/cohort_mutations_analysis/main"

workflow GENOME_INTERPRETER {
    take:
    cnaqc_out   // tuple val(meta), path(file)
    tinc_out    // tuple val(meta), path(file)
    join_cnaqc_out
    tmb_rds

    main:
    ch_versions = Channel.empty()
    summary_table_rds = null
    summary_plot_rds = null
    summary_report_pdf = null
    oncoprint = null
   
    
    // cnaqc_qc_rds = cnaqc_out
    //     .filter { meta, file ->
    //         file.name.endsWith('_qc.rds')
    //     }

    // tinc_fit_rds = tinc_out
    //     .filter { meta, file ->
    //         file.name.endsWith('_fit.rds')
    //     }

    // Group CNAqc files per cohort
    // Output shape: tuple(meta, [file1, file2, ...])

    cohort_cnaqc = cnaqc_out.map { meta, file ->
            meta = meta + [ id: "${meta.dataset}" ]
            [meta.subMap('dataset', 'id'), file]}
        .groupTuple()

    cohort_tinc = tinc_out.map { meta, file ->
            meta = meta + [ id: "${meta.dataset}" ]
            [meta.subMap('dataset', 'id'), file]}
        .groupTuple()

    // cohort_tinc = tinc_out
    //     .map { meta, file ->
    //         def cohort_meta = meta + [ id: meta.dataset ]
    //         tuple(cohort_meta.subMap('dataset', 'id'), file)
    //     }
    //     .groupTuple()

    // Join grouped CNAqc and TINC inputs by cohort meta
    // Output shape: tuple(meta, cnaqc_rds_files, tinc_rds_files)

    // commenting it out temporarly 
    cohort_qc_input = cohort_cnaqc.join(cohort_tinc)

    COHORT_QC(cohort_qc_input)

    ch_versions = ch_versions.mix(COHORT_QC.out.versions)
    summary_table_rds  = COHORT_QC.out.summary_table_rds
    summary_plot_rds  = COHORT_QC.out.summary_plot_rds
    summary_report_pdf = COHORT_QC.out.summary_report_pdf


    join_cnaqc_out = join_cnaqc_out.map{ meta, rds, samples ->
        def patient = meta.patient
        meta = meta + [id: "${meta.dataset}" ]
        [meta.subMap('dataset', 'id'), rds, patient]}
        .groupTuple()

    // join_cnaqc_out = join_cnaqc_out.map { meta, rds ->
    //     def patient = meta.patient
    //     // meta = (meta + [id: "${meta.dataset}"]).subMap(['dataset', 'id'])
    //     // [meta, rds, patient]
    //     def newMeta = [dataset: "${meta.dataset}", id: "${meta.dataset}"]
    //     meta = newMeta
    //     [meta, rds, patient]
    // }.groupTuple()

    tmb_rds = tmb_rds.map { meta, rds ->
        def patient = meta.patient
        meta = meta + [ id: "${meta.dataset}" ]
        [meta.subMap('dataset', 'id'), rds, patient]}
    .groupTuple()

    // join_cnaqc_out.view()
    // tmb_rds.view()

    oncoprint_input = join_cnaqc_out.join(tmb_rds)
    oncoprint_input.view()

    COHORT_MUTATIONS(oncoprint_input)
    ch_versions = ch_versions.mix(COHORT_MUTATIONS.out.versions)

    oncoprint = COHORT_MUTATIONS.out.cohort_oncoprint
    // summary_table_rds  = COHORT_MUTATIONS.out.summary_table_rds
    // summary_plot_rds  = COHORT_MUTATIONS.out.summary_plot_rds
    // summary_report_pdf = COHORT_MUTATIONS.out.summary_report_pdf


    emit:
    summary_table_rds
    summary_plot_rds
    summary_report_pdf
    oncoprint
    versions = ch_versions

}
