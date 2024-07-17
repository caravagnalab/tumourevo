process ANNOTATE_DRIVER {
    container = 'docker://lvaleriani/cnaqc:dev1'

    input:

    tuple val(meta), path(snv_RDS)

    output:

    tuple val(meta), path("*.rds"), emit: rds

    script:
    def args                                = task.ext.args    
    def prefix                              = task.ext.prefix                                       ?: "${meta.id}"  

    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(readr)

    data = readRDS("$snv_RDS")
    SNV = data[["$meta.tumour_sample"]]
    SNV = SNV\$mutations

    drivers_table = readr::read_tsv(file = "$params.drivers_table") 
    
    if("$meta.cancer_type" == 'PANCANCER'){
      drivers_table = drivers_table %>% 
        dplyr::group_by(SYMBOL) %>% 
        dplyr::reframe(CGC_CANCER_GENE = any(CGC_CANCER_GENE), dplyr::across(dplyr::everything())) %>% 
        dplyr::filter(CGC_CANCER_GENE) %>% 
        dplyr::mutate(CANCER_TYPE = 'PANCANCER')
    } 

    drivers_table = drivers_table %>% 
        dplyr::select(SYMBOL, CANCER_TYPE, CGC_CANCER_GENE) %>% 
        unique()


    x = SNV %>% 
      dplyr::mutate(CANCER_TYPE = "$meta.cancer_type") %>%
      dplyr::left_join(
        drivers_table,
        by = c('SYMBOL', 'CANCER_TYPE')
      ) %>% 
      tidyr::separate(HGVSp, ':', into = c('s1', 's2'), remove=F) %>%
      dplyr::mutate(
          is_driver = (CGC_CANCER_GENE & IMPACT %in% c('MODERATE', 'HIGH')),
          driver_label = paste(SYMBOL, s2)
      )

    new_data = list()
    new_data[["$meta.tumour_sample"]]\$mutations = x
    new_data[["$meta.tumour_sample"]]\$sample = data[["$meta.tumour_sample"]]\$sample
    saveRDS(object = new_data, file = paste0("$prefix", ".rds"))
    """
}
