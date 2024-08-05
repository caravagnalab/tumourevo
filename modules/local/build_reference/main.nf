//
// BUILD_REFERENCE PROCESS
//
// Build reference doc: 
//   http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/buildref.html
//
// We need 2 additional files:
//  - the path to a tab-delimited table of transcripts
//  - the path to a fasta file for the reference genome of interest
//
// For now I am using params for input files.
// Testing with files: BioMart_human_GRCh37_chr3_segment.txt and chr3_segment.fa

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
