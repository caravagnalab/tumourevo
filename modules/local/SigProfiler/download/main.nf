process DOWNLOAD_GENOME_SIGPROFILER {
    
    container = 'docker://katiad/sigprofiler:latest'

    input:
      tuple val(reference_genome) // reference_genome : download_genome_sigprofiler_reference_genome -> for example: GRCh37

    output:
      path("*"), emit: genome 
     
      
    script:
    
      //def sigprofiler_genome_path           = args!='' && args.sigprofiler_genome_path      ? "$args.sigprofiler_genome_path" : ""
      def sigprofiler_genome_path = "signature_deconvolution/Sigprofiler/genome/"

      """
      if [ ! -d "$sigprofiler_genome_path" ]; then

          mkdir -p $sigprofiler_genome_path
      else
          echo "Directory $sigprofiler_genome_path already exists."
      fi

      SigProfilerMatrixGenerator install $reference_genome -v $sigprofiler_genome_path

      """
}
                                                                                                       
