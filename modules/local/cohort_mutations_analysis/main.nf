process COHORT_MUTATIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/1099ace76d5a6d0f47a9da4a8b136906c5423ef83b591868b94a694703eeae86/data':
        'community.wave.seqera.io/library/bioconductor-complexheatmap_r-cnaqc_r-dplyr_r-ggplot2_r-tidyr:2469333aec77238e' }"
        
    input:
    tuple val(meta), path(join_cnaqc), val(cnaqc_patients), path(tmb_rds), val(tmb_patients)

    output:
    tuple val(meta), path("*_oncoprint.pdf"), emit: cohort_oncoprint
    // tuple val(meta), path("*_tmb.rds"), emit: tmb_rds 
    // tuple val(meta), path("*_vaf_chr_plot.rds"), path('*_chr_mut.rds'), path('*_consequence_mut_plot.rds'), path('*_mut_type_plot.rds'), path('*_driver_oncoprint.rds'), emit: data_stats_rds
    // tuple val(meta), path('*_mutations_report.pdf'), emit: mutation_report_pdf
    // tuple val(meta), path('*_driver_oncoprint.pdf'), emit: driver_oncoprint_plot_pdf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    template "main_script.R"

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${prefix}_oncoprint.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        complexheatmap: \$(Rscript -e "library(ComplexHeatmap); cat(as.character(packageVersion('ComplexHeatmap')))")
        patchwork: \$(Rscript -e "library(patchwork); cat(as.character(packageVersion('patchwork')))")
        cnaqc: \$(Rscript -e "library(CNAqc); cat(as.character(packageVersion('CNAqc')))")
    END_VERSIONS
    """
}
