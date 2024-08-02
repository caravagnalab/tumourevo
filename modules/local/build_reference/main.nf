process BUILD_REFERENCE {
  tag "$meta.id"
  container='file:///fast/cdslab/ebusca00/singularity/cdslab.sif'

  input:
    
    tuple val(meta), path(cds), path(genome)
  
  output:

    tuple val(meta), path("reference.rda"), emit: dnds_reference

  script:

    def args                                = task.ext.args    
    def prefix                              = task.ext.prefix                                       ?: "${meta.id}"  

    """
    #!/usr/bin/env Rscript

    library(dndscv)
    
    buildref(
      cdsfile="$cds", 
      genomefile="$genome", 
      outfile = "reference.rda", 
      excludechrs="MT"
    )
    """
}
