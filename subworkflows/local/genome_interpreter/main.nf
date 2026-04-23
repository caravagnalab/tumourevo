//
// GENOME INTERPRETER WORKFLOW
//

include { COHORT_QC } from "../../../modules/local/cohort_qc/main"
include { COHORT_MUTATIONS } from "../../../modules/local/cohort_mutations_analysis/main"
include { SUBCLONAL_INTERPRETATION } from "../../../modules/local/subclonal_interpretation/main"


workflow GENOME_INTERPRETER {
    take:
    cnaqc_out   // tuple val(meta), path(file)
    tinc_out    // tuple val(meta), path(file)
    join_cnaqc_out
    tmb_rds
    table_pyclone
    table_mobster
    table_viber
    assign_pyclone
    assign_mobster
    assign_viber

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

    // Prepare inputs for subclonal interpretation
    mutation_tables_ch = table_mobster.map {meta, table ->
        meta = meta + [id: "${meta.dataset}_${meta.patient}"]
        [meta.subMap('dataset', 'patient', 'id'), table]
    }.groupTuple(by: 0)
    .join(table_pyclone)
    .join(table_viber)
    .map { tuple ->
        def meta = tuple[0]
        def mobster_files = tuple[1]
        def pyclone_file = tuple[2]
        def viber_file = tuple[3]
        [meta, mobster_files + [pyclone_file, viber_file]]
    }

    results_sigprofiler_ch = assign_pyclone
        .join(assign_viber)
        .map { tuple ->
            def meta = tuple[0]
            def sigprofiler_files = tuple[1..-1]
            [meta, sigprofiler_files]
        }

    subclonal_input = mutation_tables_ch
        .combine(results_sigprofiler_ch, by: 0)

    SUBCLONAL_INTERPRETATION(subclonal_input)
    // summary_subclonal_pdf = SUBCLONAL_INTERPRETATION.out.report_pdf

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
    // oncoprint_input.view()

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
    // summary_subclonal_pdf
    versions = ch_versions

}
