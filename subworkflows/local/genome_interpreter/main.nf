//
// GENOME INTERPRETER WORKFLOW
//

include { COHORT_QC_REPORT } from "../../../modules/local/cohort_qc_report/main"


workflow GENOME_INTERPRETER {
    take:
    cnaqc_out   // tuple val(meta), path(file)
    tinc_out    // tuple val(meta), path(file)

    main:
    ch_versions = Channel.empty()
    summary_table_rds = null
    summary_plot_rds = null
    summary_report_pdf = null
   
    
    cnaqc_qc_rds = cnaqc_out
        .filter { meta, file ->
            file.name.endsWith('_qc.rds')
        }

    tinc_fit_rds = tinc_out
        .filter { meta, file ->
            file.name.endsWith('_fit.rds')
        }

    // Group CNAqc files per cohort
    // Output shape: tuple(meta, [file1, file2, ...])

    cohort_cnaqc = cnaqc_qc_rds
        .map { meta, file ->
            def cohort_meta = meta + [ id: meta.dataset ]
            tuple(cohort_meta.subMap('dataset', 'id'), file)
        }
        .groupTuple()

    cohort_tinc = tinc_fit_rds
        .map { meta, file ->
            def cohort_meta = meta + [ id: meta.dataset ]
            tuple(cohort_meta.subMap('dataset', 'id'), file)
        }
        .groupTuple()

    // Join grouped CNAqc and TINC inputs by cohort meta
    // Output shape: tuple(meta, cnaqc_rds_files, tinc_rds_files)

    cohort_qc_input = cohort_cnaqc.join(cohort_tinc)

    COHORT_QC_REPORT(cohort_qc_input)
    
    ch_versions = ch_versions.mix(COHORT_QC_REPORT.out.versions)
    summary_table_rds  = COHORT_QC_REPORT.out.summary_table_rds
    summary_plot_rds  = COHORT_QC_REPORT.out.summary_plot_rds
    summary_report_pdf = COHORT_QC_REPORT.out.summary_report_pdf

    emit:
    summary_table_rds
    summary_plot_rds
    summary_report_pdf
    versions = ch_versions

}
