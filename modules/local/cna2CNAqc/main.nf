process CNA_PROCESSING {
    tag "$meta.id"

    container "docker://lvaleriani/cnaqc:dev1"

    input:
      tuple val(meta), path(cnaPath)

    output:
      tuple val(meta), path("*.rds"), emit: rds

    when:
      task.ext.when == null || task.ext.when

    script:
      def args = task.ext.args ?: ''
      def prefix = task.ext.prefix ?: "${meta.id}"
      
      """
      #!/usr/bin/env Rscript 
      source(paste0("$moduleDir", '/parser_CNA.R'))

      if ("$meta.cna_caller" == 'sequenza'){
        CNA = parse_Sequenza("$meta.tumour_sample", "$cnaPath")

      } else if ("$meta.cna_caller" == 'ASCAT'){
        CNA = parse_ASCAT("$meta.tumour_sample", "$cnaPath")
        
      } else {
        stop('Copy Number Caller not supported.')
      }

      saveRDS(object = CNA, file = paste0(${prefix}, ".rds"))
      """
}
