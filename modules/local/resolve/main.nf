process RESOLVE {
    tag "$meta.id"
    label "process_long"
    
    //container 

    input:
        tuple val(meta), path(tsv_join,  stageAs: '*.tsv')

    output:
       
        path "versions.yml",                                    emit: versions

    script:
    template "main_script.R"
}
