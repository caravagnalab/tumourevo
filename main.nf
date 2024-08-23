#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/tumourevo
========================================================================================
 nf-core/tumourevo Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/tumourevo
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TUMOUREVO } from './workflows/tumourevo'
include { samplesheetToList } from 'plugin/nf-schema'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   PARSE INPUT FILE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Initialize fasta file with meta map:
fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
input = params.input ? Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json")) : Channel.empty()
drivers_table =  params.fasta ? Channel.fromPath(params.drivers_table).map{ it -> [ it ] }.collect() : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

 //
 // WORKFLOW: Run main analysis pipeline depending on type of input
 //
 workflow NFCORE_TUMOUREVO {

    take:
    input
    fasta
    drivers_table

    main:
    
    TUMOUREVO (
        input,
        fasta,
        drivers_table
    )

    emit:
    null
}

workflow {

    main:
    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_TUMOUREVO(
        input,
        fasta,
        drivers_table
    )

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
