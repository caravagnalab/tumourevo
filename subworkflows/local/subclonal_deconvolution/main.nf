// SUBCLONAL DECONVOLUTION WORKFLOW

include { MOBSTER } from "../../../modules/nf-core/mobster/main"
include { VIBER } from "../../../modules/nf-core/viber/main"
include { PYCLONEVI } from "../../../modules/nf-core/pyclonevi/main"
include { FORMATTER } from "../../../subworkflows/local/formatter/main"
include { CTREE as CTREE_MOBSTER } from "../../../modules/nf-core/ctree/main"
include { CTREE as CTREE_PYCLONEVI } from "../../../modules/nf-core/ctree/main"
include { CTREE as CTREE_VIBER } from "../../../modules/nf-core/ctree/main"
include { PLOT_DECONVOLUTION as PLOT_DECONVOLUTION_MOBSTER } from "../../../modules/local/plot_deconvolution/main"
include { PLOT_DECONVOLUTION as PLOT_DECONVOLUTION_VIBER } from "../../../modules/local/plot_deconvolution/main"
include { PLOT_DECONVOLUTION as PLOT_DECONVOLUTION_PYCLONE } from "../../../modules/local/plot_deconvolution/main"


workflow SUBCLONAL_DECONVOLUTION {
    take:
    rds_join // tuple val(meta), path("*.rds"), val(tumour_samples), emit: rds

    main:
    ch_versions = Channel.empty()
    mobster_pdf = null
    ctree_mobster_pdf = null
    ctree_mobster_rds = null
    viber_pdf = null
    ctree_viber_pdf = null
    ctree_viber_rds = null
    pyclone_fits = null
    pyclone_best = null
    pyclone_table = null
    ctree_pyclone_pdf = null
    ctree_pyclone_rds = null
    viber_results = null
    mobster_results = null
    viber_results_heuristic = null

    if (params.tools && params.tools.split(",").contains("mobster")) {
        joinCNAqc = rds_join.transpose().map{ meta, rds, sample ->
            meta = meta + ["tumour_sample": sample, "id":"${meta.dataset}_${meta.patient}_${sample}"]
            [meta, rds]}

        MOBSTER(joinCNAqc)
        ch_versions = ch_versions.mix(MOBSTER.out.versions_mobster)

        CTREE_MOBSTER(MOBSTER.out.mobster_best_rds)
        ch_versions = ch_versions.mix(CTREE_MOBSTER.out.versions)

        mobster_results = MOBSTER.out.mobster_best_rds
        mobster_pdf = MOBSTER.out.mobster_report_pdf
        ctree_mobster_pdf = CTREE_MOBSTER.out.ctree_report_pdf
        ctree_mobster_rds = CTREE_MOBSTER.out.ctree_rds


        join_mobster = mobster_results.map { meta, rds ->
          def key = meta.subMap('dataset', 'patient') + [id: "${meta.dataset}_${meta.patient}"]
          [key, rds]
        }.groupTuple()

        input_plot = join_mobster.map { meta, fit ->
          [meta, fit, [], [], []]
        }
        PLOT_DECONVOLUTION_MOBSTER(input_plot)
    }

    if (params.tools && params.tools.split(",").contains("viber")) {
        VIBER(rds_join)
        ch_versions = ch_versions.mix(VIBER.out.versions_viber)

        CTREE_VIBER(VIBER.out.viber_rds)
        ch_versions = ch_versions.mix(CTREE_VIBER.out.versions)

        viber_results = VIBER.out.viber_rds
        viber_results_heuristic = VIBER.out.viber_heuristic_rds
        viber_pdf = VIBER.out.viber_report_pdf
        ctree_viber_pdf = CTREE_VIBER.out.ctree_report_pdf
        ctree_viber_rds = CTREE_VIBER.out.ctree_rds

        input_plot = viber_results.map { meta, fit ->
          [meta, [], fit, [], []]
        }
        PLOT_DECONVOLUTION_VIBER(input_plot)
    }

    if (params.tools && params.tools.split(",").contains("pyclone-vi")) {
        FORMATTER(rds_join, "rds")
        PYCLONEVI(FORMATTER.out.out_data)
        CTREE_PYCLONEVI(PYCLONEVI.out.ctree_input)

        ch_versions = ch_versions.mix(FORMATTER.out.versions)
        ch_versions = ch_versions.mix(PYCLONEVI.out.versions_pyclonevi)
        ch_versions = ch_versions.mix(CTREE_PYCLONEVI.out.versions)

        pyclone_fits = PYCLONEVI.out.pyclone_all_fits
        pyclone_best = PYCLONEVI.out.pyclone_best_fit
        pyclone_table = FORMATTER.out.out_data
        ctree_pyclone_pdf = CTREE_PYCLONEVI.out.ctree_report_pdf
        ctree_pyclone_rds = CTREE_PYCLONEVI.out.ctree_rds

        pyclone_table = pyclone_table.map{meta, result, sample -> [meta, result]}
        input_plot = pyclone_table.join(pyclone_best).map { meta, table, fit ->
          [meta, [], [], table, fit]
        }
        PLOT_DECONVOLUTION_PYCLONE(input_plot)
    }

    emit:
    pyclone_fits
    pyclone_best
    ctree_pyclone_pdf
    ctree_pyclone_rds
    viber_pdf
    ctree_viber_pdf
    ctree_viber_rds
    mobster_pdf
    ctree_mobster_pdf
    ctree_mobster_rds
    pyclone_table
    mobster_results
    viber_results
    viber_results_heuristic
    versions = ch_versions
}
