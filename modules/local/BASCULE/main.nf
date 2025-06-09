process BASCULE {
    tag "$meta.id"
    label "process_single"
    label "error_ignore"
    // container  // ADD

    input:
        tuple val(meta), path(tsv_join,  stageAs: '*.tsv')
    
    // output:  // ADD

}
