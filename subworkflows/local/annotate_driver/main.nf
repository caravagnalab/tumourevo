//
// DRIVER_ANNOTATION SUB-WORKFLOW
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

include { BUILD_REFERENCE } from '../../../modules/local/annotate_driver/main'
include { DNDSCV } from '../../../modules/local/annotate_driver/main'


workflow DRIVER_ANNOTATION {
    take:
        rds
        driver_list
    
    main:
        ref = BUILD_REFERENCE() // create reference for variant annotation
        DNDSCV(rds, driver_list, ref) // add dndsCV statistics and columns "known_driver" (based on driver list, default IntoGen) and "potential_driver"

    emit:
        ANNOTATE_DRIVER.out.rds

}
