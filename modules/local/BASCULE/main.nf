process BASCULE {
    tag "$meta.id"
    label "process_single"
    label "error_retry"
    // container  // ADD

    input:
        tuple val(meta), path(tsv_join,  stageAs: '*.tsv')

    output:
        tuple val(meta), path("*bascule_fit.rds"), emit: bascule_rds
        tuple val(meta), path("*_plots_all.pdf"), emit: bascule_plots_pdf
        tuple val(meta), path("*_plots_all.rds"), emit: bascule_plots_rds
        path "versions.yml", emit: versions

    script:
    template "main_script.R"

}
