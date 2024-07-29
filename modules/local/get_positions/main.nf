process GET_POSITIONS {
    tag "$meta.id"
    container = 'docker://lvaleriani/cnaqc:dev1'

    input:

    tuple val(meta), path(rds_list, stageAs: '*.rds') 

    output:

    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.bed"), emit: pos

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    res_dir = paste0("$datasetID", "/", "$patientID", "/")
    dir.create(res_dir, recursive = TRUE)
    
    samples = substr("$sampleID", 2, nchar("$sampleID")-1)
    samples = strsplit(samples, ", ")[[1]]

    positions = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(rds){
        df = readRDS(rds)
        df = df[[1]]\$mutations %>% dplyr::mutate(id = paste(chr, from, to, sep = ":")) %>% dplyr::select(chr, from, to, ref, alt, id)
    }) 

    all_positions = positions %>% dplyr::bind_rows() %>% dplyr::pull(id) %>% unlist() %>% unique()
    all = positions %>% dplyr::bind_rows() %>% dplyr::distinct() %>% dplyr::select(-id)

    dir.create(paste0("lifter/mpileup/", res_dir), recursive = T, showWarnings = F)
    write.table(file = paste0("lifter/mpileup/", res_dir, 'all_positions.bed'), all, quote = F, sep = "\t", row.names = F, col.names = T)

    missed = lapply(seq(1,length(positions)), FUN = function(s){
        all_positions[!(all_positions %in% positions[[s]]\$id)]
    }) 

    for (i in seq(1,length(missed))){
        sample = samples[[i]]
        df = dplyr::tibble(id = missed[[i]]) %>% 
                tidyr::separate(id,into = c('chr', 'from', 'to'), sep = ':') %>% 
                dplyr::filter(chr %in% c(paste0('chr', seq(1,22)), 'chrX', 'chrY'))

        write.table(file = paste0("lifter/mpileup/", res_dir, sample, "/positions_missing.bed"), df, quote = F, sep = "\t", row.names = F, col.names = F)
    }
    """
}
