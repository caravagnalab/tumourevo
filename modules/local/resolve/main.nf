process RESOLVE {
    tag "$meta.id"
    label "process_long"
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e35711a7325309ac410d78057329ec42abc40732:ba3901d0fb2269789f08df6de1e46530140dc392-0':
        'biocontainers/mulled-v2-e35711a7325309ac410d78057329ec42abc40732:ba3901d0fb2269789f08df6de1e46530140dc392-0 ' }"

    input:
        tuple val(meta), path(tsv_join,  stageAs: '*.tsv')

    output:
        tuple val(meta), path("*_mut_counts.rds")             , emit: signatures_mutCounts_rds
        tuple val(meta), path("*_fit_results.rds")            , emit: signatures_fit_results  
        path "versions.yml"                                   , emit: versions


    when:
    task.ext.when == null || task.ext.when


    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mut_counts.rds
    touch ${prefix}_fit_results.rds
  

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-resolve: \$(Rscript -e "library(RESOLVE); cat(as.character(packageVersion('RESOLVE')))")
        bioconductor-bsgenome.hsapiens.1000genomes.hs37d5: \$(Rscript -e "library(BSgenome.Hsapiens.1000genomes.hs37d5); cat(as.character(packageVersion('BSgenome.Hsapiens.1000genomes.hs37d5')))")
        bioconductor-bsgenome.hsapiens.ucsc.hg38: \$(Rscript -e "library(BSgenome.Hsapiens.UCSC.hg38); cat(as.character(packageVersion('BSgenome.Hsapiens.UCSC.hg38')))")
    END_VERSIONS
    """
}
