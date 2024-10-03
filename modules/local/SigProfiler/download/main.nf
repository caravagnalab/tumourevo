process DOWNLOAD_GENOME_SIGPROFILER {
    
    container = 'docker://katiad/sigprofiler:latest'

    input:
      val(reference_genome) // reference_genome : download_genome_sigprofiler_reference_genome -> for example: GRCh37

    output:
      path("*"), emit: genome_sigprofiler 
     
      
    script:
  

      """

      SigProfilerMatrixGenerator install $reference_genome -v .

      """
}
                                                                                                       
