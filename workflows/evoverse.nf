//include { VARIANT_ANNOTATION } from "${baseDir}/subworkflows/variant_annotation/main"
include { VCF_ANNOTATE_ENSEMBLVEP } from '../subworkflows/nf-core/vcf_annotate_ensemblvep/main'

include { FORMATTER as FORMATTER_CNA } from "${baseDir}/subworkflows/local/formatter/main"
include { FORMATTER as FORMATTER_VCF} from "${baseDir}/subworkflows/local/formatter/main"
include { LIFTER } from "${baseDir}/subworkflows/local/lifter/main"
include { DRIVER_ANNOTATION } from "${baseDir}/subworkflows/local/annotate_driver/main"
include { FORMATTER as FORMATTER_RDS} from "${baseDir}/subworkflows/local/formatter/main"
include { QC } from "${baseDir}/subworkflows/local/QC/main"
include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/local/subclonal_deconvolution/main"
include { SIGNATURE_DECONVOLUTION } from "${baseDir}/subworkflows/local/signature_deconvolution/main"
include { PLOT_REPORT_SINGLE_SAMPLE } from "${baseDir}/modules/plot_report/main"
include { PLOT_REPORT_MULTI_SAMPLE } from "${baseDir}/modules/plot_report/plot_report_multi"

workflow EVOVERSE {

  take:
  input_vcf
  input_cna
  cancer_type
  
  main:
 
  VARIANT_ANNOTATION(input_vcf)
  FORMATTER_VCF(VARIANT_ANNOTATION.out.vep, "vcf")
  FORMATTER_CNA(input_cna, "cna")

  exist_bam_val = false //placeholder
  if (params.mode == 'multisample' && exist_bam_val){
    tumor_bam = Channel.fromPath(params.samples).
      splitCsv(header: true).
       map{row ->
       tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.tumour_bam), file(row.tumour_bai))}

    LIFTER(FORMATTER_VCF.out, tumor_bam)
    annotation = DRIVER_ANNOTATION(LIFTER.out, cancer_type)

  } else {
   annotation = DRIVER_ANNOTATION(FORMATTER_VCF.out, cancer_type)
  }

  QC(FORMATTER_CNA.out, annotation)
  SUBCLONAL_DECONVOLUTION(QC.out.rds_join)
  SIGNATURE_DECONVOLUTION(QC.out.rds_join)
  
  emit:

  null
}
