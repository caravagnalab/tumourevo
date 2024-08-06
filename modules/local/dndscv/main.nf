//
// DNDSCV PROCESS
//
// dndscv module based on knitr prototype
// dndscv script (SINGLE SAMPLE ONLY)
// TODO multisamples
// 1. load input rds get mutations and sample id  
// 2. create dndscv input
// 3. Get list of drivers from file
// 4. Use custom reference

process DNDSCV {
  debug true
  tag "$meta.id"
  container "${workflow.containerEngine == 'singularity' ? 'docker://tucano/dndscv:latest' : 'tucano/dndscv:latest'}"


  input:
    
    tuple val(meta), path(snv_rds), path(driver_list), path(reference)
  
  output:

    tuple val(meta), path("*_dnds.rds"), emit: dnds_rds
  
  script:

    def args                                = task.ext.args ?: ""
    def prefix                              = task.ext.prefix ?: "${meta.id}"  
    
    """
    dndscv_runner.R \\
      -i ${snv_rds} \\
      -s ${meta.tumour_sample} \\
      -d ${driver_list} \\
      -r ${reference} \\
      -o "${prefix}_dnds.rds" \\
      --qmis_cv ${params.dndscv_qmis_cv} \\
      --qtrunc_cv ${params.dndscv_qtrunc_cv} \\
      --qallsubs_cv ${params.dndscv_qallsubs_cv} \\
      ${args}
    """
}
