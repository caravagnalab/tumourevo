//
// GENOME INTERPRETER WORKFLOW
//

include { COHORT_QC } from "../../../modules/local/cohort_qc/main"
// include { COHORT_MUTATIONS } from "../../../modules/local/cohort_mutations_analysis/main"
include { SUBCLONAL_INTERPRETATION } from "../../../modules/local/subclonal_interpretation/main"
include { COHORT_SIGNATURES } from "../../../modules/local/cohort_signatures/main"

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
    sigprofiler_out
    sparsesignature_assign_cosmic
    

    main:
    ch_versions = Channel.empty()
    summary_table_rds = null
    summary_plot_rds = null
    summary_report_pdf = null
    summary_subclonal_pdf = null
    summary_subclonal_rds = null
    // oncoprint = null
    cohort_signatures_pdf = null
    cohort_signatures_rds = null

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

    // Prepare inputs for subclonal interpretation
    if (params.tools && params.tools.split(",").contains("mobster") && params.tools.split(",").contains("sigprofiler")) {
        if (params.tools && params.tools.split(",").contains("viber") && params.tools.split(",").contains("pyclone-vi")) {
            mutation_tables_ch = table_mobster.map { meta, table ->
                meta = meta + [id: "${meta.dataset}_${meta.patient}"]
                [meta.subMap('dataset', 'patient', 'id'), table]
            }.groupTuple(by: 0)
            .join(table_pyclone)
            .join(table_viber)
            .map { meta, mobster_files, pyclone_file, viber_file ->
                [meta, mobster_files + [pyclone_file, viber_file]]
            }
        } else if (params.tools && params.tools.split(",").contains("viber") && !params.tools.split(",").contains("pyclone-vi")) {
            mutation_tables_ch = table_mobster.map { meta, table ->
                meta = meta + [id: "${meta.dataset}_${meta.patient}"]
                [meta.subMap('dataset', 'patient', 'id'), table]
            }.groupTuple(by: 0)
            .join(table_viber)
            .map { meta, mobster_files, viber_file ->
                [meta, mobster_files + [viber_file]]
            }
        } else if (params.tools && params.tools.split(",").contains("pyclone-vi") && !params.tools.split(",").contains("viber")) {
            mutation_tables_ch = table_mobster.map { meta, table ->
                meta = meta + [id: "${meta.dataset}_${meta.patient}"]
                [meta.subMap('dataset', 'patient', 'id'), table]
            }.groupTuple(by: 0)
            .join(table_pyclone)
            .map { meta, mobster_files, pyclone_file ->
                [meta, mobster_files + [pyclone_file]]
            }
        }

        results_sigprofiler_ch = assign_pyclone
            .join(assign_viber)
            .map { meta, file1, file2 ->
                [meta, [file1, file2]]
            }

        subclonal_input = mutation_tables_ch
            .combine(results_sigprofiler_ch, by: 0)

        SUBCLONAL_INTERPRETATION(subclonal_input)
        summary_subclonal_pdf = SUBCLONAL_INTERPRETATION.out.report_pdf
        summary_subclonal_rds = SUBCLONAL_INTERPRETATION.out.rds
    }

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

    // disabling momentarly the cohort mutations analysis and visualization
    // COHORT_MUTATIONS(oncoprint_input)
    // ch_versions = ch_versions.mix(COHORT_MUTATIONS.out.versions)

    // oncoprint = COHORT_MUTATIONS.out.cohort_oncoprint
    // summary_table_rds  = COHORT_MUTATIONS.out.summary_table_rds
    // summary_plot_rds  = COHORT_MUTATIONS.out.summary_plot_rds
    // summary_report_pdf = COHORT_MUTATIONS.out.summary_report_pdf

    if (params.tools && params.tools.split(",").contains("sigprofiler") && !params.tools.split(",").contains("sparsesignatures")) {
        cohort_sigprofiler = sigprofiler_out.map { meta, file ->
            meta = meta + [ id: "${meta.dataset}" ]
            [meta.subMap('dataset', 'id'), file]}
        cohort_signatures_input = cohort_sigprofiler
    } else if (params.tools && params.tools.split(",").contains("sparsesignatures") && !params.tools.split(",").contains("sigprofiler")) {
        cohort_sparsesig = sparsesignature_assign_cosmic.map { meta, file ->
            meta = meta + [ id: "${meta.dataset}" ]
            [meta.subMap('dataset', 'id'), file]}
        cohort_signatures_input = cohort_sparsesig
    } else {
        cohort_sigprofiler = sigprofiler_out.map { meta, file ->
            meta = meta + [ id: "${meta.dataset}" ]
            [meta.subMap('dataset', 'id'), file]}
        cohort_sparsesig = sparsesignature_assign_cosmic.map { meta, file ->
            meta = meta + [ id: "${meta.dataset}" ]
            [meta.subMap('dataset', 'id'), file]}
        cohort_signatures_input  = cohort_sigprofiler.join(cohort_sparsesig, remainder: true).map { meta, file1, file2 ->
            [meta, [file1, file2]]}.map { meta, file -> [meta, file.flatten()]}
    }
    
    cohort_signatures_input.view()
    COHORT_SIGNATURES(cohort_signatures_input)
    cohort_signatures_pdf = COHORT_SIGNATURES.out.report_cohort_signatures
    cohort_signatures_rds = COHORT_SIGNATURES.out.rds_cohort_signatures
    ch_versions = ch_versions.mix(COHORT_SIGNATURES.out.versions)
    
    emit:
    summary_table_rds
    summary_plot_rds
    summary_report_pdf
    // oncoprint
    summary_subclonal_pdf
    summary_subclonal_rds
    cohort_signatures_pdf
    cohort_signatures_rds
    versions = ch_versions

}
