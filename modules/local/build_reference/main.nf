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
// 
// I do not need a meta here
//
process BUILD_REFERENCE {

  tag "BUILD_REFERENCE"
  container "${workflow.containerEngine == 'singularity' ? 'docker://tucano/dndscv:latest' : 'tucano/dndscv:latest'}"

  input:
    path(cds)
    path(genome)
  
  output:
    path("reference.rda"), emit: dnds_reference

  script:
    // no nedd for arg or prefix here    
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
