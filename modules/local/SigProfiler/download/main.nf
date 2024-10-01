process DOWNLOAD_GENOME_SIGPROFILER {
    //tag "$meta.id"
    container = 'docker://katiad/sigprofiler:latest'

    input:
      tuple val(reference_genome) // reference_genome : download_genome_sigprofiler_reference_genome -> for example: GRCh37

    output:
      path("*"), emit: genome 
     
      
    script:
    
      //def args                              = task.ext.args                                 ?: ''
      //def prefix                            = task.ext.prefix                               ?: "${meta.id}"
      //def reference_genome                  = args!='' && args.reference_genome             ? "$args.reference_genome" : ""
      //def sigprofiler_genome_path           = args!='' && args.sigprofiler_genome_path      ? "$args.sigprofiler_genome_path" : ""
      //def sigprofiler_genome_path = "$meta.datasetID/signature_deconvolution/Sigprofiler/genome/"

      """
      SigProfilerMatrixGenerator install $reference_genome -v .
      """
}
                                                                                                       
