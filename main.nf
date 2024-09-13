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

include { ANNOTATION_CACHE_INITIALISATION  } from './subworkflows/local/annotation_cache_initialisation'
include { DOWNLOAD_CACHE_VEP } from './subworkflows/local/download_cache_vep'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   PARSE INPUT FILE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

input = params.input ? Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json")) : Channel.empty()

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

    main:
    fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
    drivers_table =  params.fasta ? Channel.fromPath(params.drivers_table).map{ it -> [ it ] }.collect() : Channel.empty()

    if (params.download_cache_vep) {
        ensemblvep_info = Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
        DOWNLOAD_CACHE_VEP(ensemblvep_info)
        vep_cache = DOWNLOAD_CACHE_VEP.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }
        
    } else {
        ANNOTATION_CACHE_INITIALISATION(
            params.vep_cache,
            params.vep_species,
            params.vep_cache_version,
            params.vep_genome,
            params.vep_custom_args)

        vep_cache  = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache
    }
    
    TUMOUREVO (
        input,
        fasta,
        drivers_table,
        vep_cache
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
        input
    )

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
