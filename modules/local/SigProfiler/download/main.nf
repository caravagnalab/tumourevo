process DOWNLOAD_SIGPROF {
    
    publishDir params.publish_dir, mode: 'copy'

    input:
       //tuple val(datasetID), val(patientID), val(sampleID)
       tuple val(meta)

    output:
       tuple val(datasetID), path("signature_deconvolution/Sigprofiler/"), emit: genome 
     
      
    script:
    
      def args                              = task.ext.args                                 ?: ''
      def reference_genome                  = args!='' && args.reference_genome             ? "$args.reference_genome" : ""
      


      """

      mkdir -p signature_deconvolution/Sigprofiler/genome
  
      SigProfilerMatrixGenerator install $reference_genome -v signature_deconvolution/Sigprofiler/genome/
      
      
      """
}
                                                                                                               
