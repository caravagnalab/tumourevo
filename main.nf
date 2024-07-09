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

//include { TUMOUREVO } from './workflows/tumourevo'
include { samplesheetToList } from 'plugin/nf-schema'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   PARSE INPUT FILE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Initialize fasta file with meta map:
fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
input = params.input ? Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json")) : Channel.empty()
input.view()

//input_vcf = Channel.fromPath(params.input).
//    splitCsv(header: true).
//    map {
//      row ->
//      meta = row.subMap('dataset', 'patient', 'tumour_sample', 'normal_sample', 'cancer_type', 'cnv_caller')
//      [meta, [
//          file(row.vcf),
//          file(row.vcf_tbi),
//          file(row.tumour_bam)
//          file(row.tumour_bai) 
//      ]]
//    }

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

    main:

    TUMOUREVO (
        input,
        fasta
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
        fasta
    )

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
