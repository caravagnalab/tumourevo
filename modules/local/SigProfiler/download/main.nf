process DOWNLOAD_GENOME_SIGPROFILER {
    tag "$meta.id"
    container = 'docker://katiad/sigprofiler:latest'
   

    input:
       //tuple val(datasetID), val(patientID), val(sampleID)
       tuple val(meta)

    output:
       tuple val(datasetID), path("signature_deconvolution/Sigprofiler/"), emit: genome 
     
      
    script:
    
      def args                              = task.ext.args                                 ?: ''
      def reference_genome                  = args!='' && args.reference_genome             ? "$args.reference_genome" : ""
      def sigprofiler_genome_path = "signature_deconvolution/Sigprofiler/genome"


      """
      
      if [ ! -d "$sigprofiler_genome_path" ]; then

          mkdir -p $sigprofiler_genome_path
      else
          echo "Directory $sigprofiler_genome_path already exists."
      fi
      
  
      SigProfilerMatrixGenerator install $reference_genome -v $sigprofiler_genome_path
      

      """

}
                                                                                                               
