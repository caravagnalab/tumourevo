process COHORT_QC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://community.wave.seqera.io/library/bioconductor-complexheatmap_r-circlize_r-gridextra_r-patchwork_pruned:2e84eed749dc5257' :
        'community.wave.seqera.io/library/bioconductor-complexheatmap_r-circlize_r-gridextra_r-patchwork_pruned:2e84eed749dc5257' }"

    input:
    tuple val(meta), path(cnaqc_rds_files), path(tinc_rds_files)

    output:
    tuple val(meta), path("*.qc_summary.rds"),   emit: summary_table_rds
    tuple val(meta), path("*.cna_segments.rds"), emit: summary_cna_segments_rds
    tuple val(meta), path("*.qc_plot.rds"),      emit: summary_plot_rds
    tuple val(meta), path("*.qc_report.pdf"),    emit: summary_report_pdf
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.qc_summary.rds
    touch ${prefix}.qc_plot.rds
    touch ${prefix}.qc_report.pdf
    touch ${prefix}.cna_segments.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        ComplexHeatmap: \$(Rscript -e "cat(as.character(packageVersion('ComplexHeatmap')))")
        circlize: \$(Rscript -e "cat(as.character(packageVersion('circlize')))")
        scales: \$(Rscript -e "cat(as.character(packageVersion('scales')))")
        patchwork: \$(Rscript -e "cat(as.character(packageVersion('patchwork')))")
    END_VERSIONS
    """
}


