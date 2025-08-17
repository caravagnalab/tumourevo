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
        
        path "versions.yml",                                    emit: versions

    script:
    template "main_script.R"
}
