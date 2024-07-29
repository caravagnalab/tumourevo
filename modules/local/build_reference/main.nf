process BUILD_REFERENCE {
  tag "$meta.id"
  container = '<container_path>'

  input:
    
    tuple val(meta)
  
  output:

    tuple val(meta), path(REFERENCE), emit: dnds_rds

  script:

    def args                                = task.ext.args    
    def prefix                              = task.ext.prefix                                       ?: "${meta.id}"  

    """
    #!/usr/bin/env Rscript


    """
}
