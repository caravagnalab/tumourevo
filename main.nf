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

input_vcf = Channel.fromPath(params.input).
    splitCsv(header: true).
    map {
      row ->
      tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.vcf), file(row.vcf_tbi))
    }

cancer_type = Channel.fromPath(params.input).
    splitCsv(header: true).
    map {
      row ->
      row.cancer_type.toString()
    }

<<<<<<< HEAD
// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-evoverse-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/evoverse Workflow Summary'
    section_href: 'https://github.com/nf-core/evoverse'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

=======
input_cna = Channel.fromPath(params.input).
  splitCsv(header: true).
  map {
    row ->
   tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.cnv_res), row.cnv_caller.toString())
  }
>>>>>>> 00e5e3d1cdf7d5bd5b950651c6610929c1d8822f

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

    main:

    EVOVERSE (
        input_vcf,
        input_cna
        cancer_type
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
    )

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
