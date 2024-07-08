#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/evoverse
========================================================================================
 nf-core/evoverse Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/evoverse
----------------------------------------------------------------------------------------
*/

*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EVOVERSE } from '${baseDir}/workflows/evoverse'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   PARSE INPUT FILE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Initialize fasta file with meta map:
fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

input_vcf = Channel.fromPath(params.input).
    splitCsv(header: true).
    map {
      row ->
      meta = row.subMap('dataset', 'patient', 'tumour_sample', 'normal_sample', 'cancer_type', 'cnv_caller')
      [meta, [
          file(row.vcf),
          file(row.vcf_tbi),
          file(row.tumour_bam)
          file(row.tumour_bai) 
      ]]
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

 //
 // WORKFLOW: Run main analysis pipeline depending on type of input
 //
 workflow NFCORE_EVOVERSE {

    take:
    input_vcf
    input_cna
    cancer_type
    fasta

    main:

    EVOVERSE (
        input_vcf,
        input_cna
        cancer_type
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
    NFCORE_EVOVERSE(
        input_vcf,
        input_cna
        cancer_type
        fasta
    )

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
