process BUILD_REFERENCE {
  tag "$meta.id"
  container = 'docker://lvaleriani/cnaqc:dev1'

  input:
    
    tuple val(meta), path(cds), path(genome)
  
  output:

    tuple val(meta), path("reference.rda"), emit: dnds_rds

  script:

    def args                                = task.ext.args    
    def prefix                              = task.ext.prefix                                       ?: "${meta.id}"  

    """
    #!/usr/bin/env Rscript

    library(dndscv)
    print(dir())
    buildref(
      cdsfile="$cds", 
      genomefile="$genome", 
      outfile = "reference.rda", 
      excludechrs="MT"
    )
    """
}
