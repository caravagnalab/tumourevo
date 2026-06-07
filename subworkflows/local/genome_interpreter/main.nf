//
// GENOME INTERPRETER WORKFLOW
//

include { COHORT_QC } from "../../../modules/local/cohort_qc/main"
include { COHORT_MUTATIONS } from "../../../modules/local/cohort_mutations_analysis/main"
include { SUBCLONAL_INTERPRETATION } from "../../../modules/local/subclonal_interpretation/main"
include { COHORT_SIGNATURES } from "../../../modules/local/cohort_signatures/main"
include { PLOT_CLONE_TREE } from "../../../modules/local/plot_clone_tree/main"

workflow GENOME_INTERPRETER {
    take:
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
    ctree_viber
    ctree_pyclone
    

    main:
    ch_versions = Channel.empty()
    summary_table_rds = null
    summary_cna_segments_rds = null
    summary_plot_rds = null
    summary_report_pdf = null
    report_score = null
    report_signature = null
    oncoprint = null
    cohort_signatures_pdf = null
    cohort_signatures_rds = null


    cohort_cnaqc = join_cnaqc_out.map { meta, file, samples ->
            meta = meta + [ id: "${meta.dataset}" ]
            [meta.subMap('dataset', 'id'), file]
    }
    .groupTuple()

    cohort_tinc = tinc_out.map { meta, file ->
            meta = meta + [ id: "${meta.dataset}" ]
            [meta.subMap('dataset', 'id'), file]}
        .groupTuple()

    cohort_qc_input = cohort_cnaqc.join(cohort_tinc)

    COHORT_QC(cohort_qc_input)
    ch_versions = ch_versions.mix(COHORT_QC.out.versions)
    summary_table_rds  = COHORT_QC.out.summary_table_rds
    summary_plot_rds  = COHORT_QC.out.summary_plot_rds
    summary_report_pdf = COHORT_QC.out.summary_report_pdf
    summary_cna_segments_rds = COHORT_QC.out.summary_cna_segments_rds


    if (params.tools && params.tools.split(",").contains("mobster") && params.tools.split(",").contains("sigprofiler")) {
      def requested_tools = params.tools.split(",")

      ch_mobster = requested_tools.contains("mobster")
          ? table_mobster.map { meta, table ->
              def key = meta.subMap('dataset', 'patient') + [id: "${meta.dataset}_${meta.patient}"]
              [key, table]
            }.groupTuple(by: 0)
          : channel.empty()

      ch_pyclone = requested_tools.contains("pyclone-vi")
          ? table_pyclone.map { meta, table ->
              [meta, [table]]
            }
          : channel.empty()

      ch_viber = requested_tools.contains("viber")
          ? table_viber.map { meta, table ->
              [meta, [table]]
            }
          : channel.empty()

      mutation_tables_ch = ch_mobster              
          .map { meta, files -> [meta, files] }           // [key, [file1, file2, ...]]
          .mix(ch_pyclone)                                // [key, [file]]
          .mix(ch_viber)                                  // [key, [file]]
          .groupTuple(by: 0)                              // [key, [[files...], [file], [file]]]
          .map { meta, file_lists ->
              [meta, file_lists.flatten()]                // [key, [all files flat]]
          }

      ch_assign_pyclone = requested_tools.contains("pyclone-vi")
          ? assign_pyclone.map { meta, file -> [meta, [file]] }
          : channel.empty()

      ch_assign_viber = requested_tools.contains("viber")
          ? assign_viber.map { meta, file -> [meta, [file]] }
          : channel.empty()

      results_sigprofiler_ch = ch_assign_pyclone
          .mix(ch_assign_viber)
          .groupTuple(by: 0)
          .map { meta, file_lists ->
              [meta, file_lists.flatten()]
          }

      subclonal_input = mutation_tables_ch
          .combine(results_sigprofiler_ch, by: 0)

      SUBCLONAL_INTERPRETATION(subclonal_input)
      report_score = SUBCLONAL_INTERPRETATION.out.report_score
      report_signature = SUBCLONAL_INTERPRETATION.out.report_signature

      input_clone_tree = SUBCLONAL_INTERPRETATION.out.rds_score
        .join(SUBCLONAL_INTERPRETATION.out.rds_signature)
        .join(ctree_viber, remainder: true)
        .join(ctree_pyclone)
        .map { tuple ->
            def meta = tuple[0]
            def rds_score = tuple[1]
            def rds_sig = tuple[2]
            def viber = tuple[3] ?: []
            def pyclone = tuple[4]
            [meta, rds_score, rds_sig, viber, pyclone]
        }
      PLOT_CLONE_TREE(input_clone_tree)
  }

    join_cnaqc_out = join_cnaqc_out.map{ meta, rds, samples ->
        def patient = meta.patient
        meta = meta + [id: "${meta.dataset}" ]
        [meta.subMap('dataset', 'id'), rds, patient]}
        .groupTuple()


    tmb_rds = tmb_rds.map { meta, rds ->
        def patient = meta.patient
        meta = meta + [ id: "${meta.dataset}" ]
        [meta.subMap('dataset', 'id'), rds, patient]}
    .groupTuple()

    oncoprint_input = join_cnaqc_out.join(tmb_rds)
   

    COHORT_MUTATIONS(oncoprint_input)
    ch_versions = ch_versions.mix(COHORT_MUTATIONS.out.versions)

    oncoprint = COHORT_MUTATIONS.out.cohort_oncoprint
    

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
    } else if (params.tools && params.tools.split(",").contains("sigprofiler") && params.tools.split(",").contains("sparsesignatures")) {

        cohort_sigprofiler = sigprofiler_out.map { meta, file ->
            meta = meta + [ id: "${meta.dataset}" ]
            [meta.subMap('dataset', 'id'), file]
        }

        cohort_sparsesig = sparsesignature_assign_cosmic.map { meta, file ->
            meta = meta + [ id: "${meta.dataset}" ]
            [meta.subMap('dataset', 'id'), file]
        }

        cohort_signatures_input = cohort_sigprofiler
            .join(cohort_sparsesig, remainder: true)
            .map { meta, file1, file2 ->
                [meta, [file1, file2].findAll { it != null }.flatten()]
            }

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
    summary_cna_segments_rds
    oncoprint
    report_score
    report_signature
    cohort_signatures_pdf
    cohort_signatures_rds
    versions = ch_versions

}
