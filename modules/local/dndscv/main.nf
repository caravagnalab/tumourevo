process DNDSCV {
  tag "$meta.id"
  container = 'quay.io/cellgeni/cellgeni-jupyter'

  input:
    
    tuple val(meta), path(snv_rds), path(driver_list), path(reference)
  
  output:

    tuple val(meta), path("*_dnds.rds"), emit: dnds_rds

  script:

    def args                                = task.ext.args    
    def prefix                              = task.ext.prefix                                       ?: "${meta.id}"  

    """
    #!/usr/bin/env Rscript

    library(stringr)
    library(dndscv)

    input <- readRDS($snv_rds)
    sample_id <- input["${meta.tumour_sample}"]["sample"]

    """
}
