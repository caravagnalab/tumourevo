include { VCF_ANNOTATE_ENSEMBLVEP } from '../subworkflows/nf-core/vcf_annotate_ensemblvep/main'
include { FORMATTER as FORMATTER_CNA } from "${baseDir}/subworkflows/local/formatter/main"
include { FORMATTER as FORMATTER_VCF} from "${baseDir}/subworkflows/local/formatter/main"

include { LIFTER } from "${baseDir}/subworkflows/local/lifter/main"
include { DRIVER_ANNOTATION } from "${baseDir}/subworkflows/local/annotate_driver/main"
include { FORMATTER as FORMATTER_RDS} from "${baseDir}/subworkflows/local/formatter/main"
include { QC } from "${baseDir}/subworkflows/local/QC/main"
include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/local/subclonal_deconvolution/main"
include { SIGNATURE_DECONVOLUTION } from "${baseDir}/subworkflows/local/signature_deconvolution/main"
//include { PLOT_REPORT_SINGLE_SAMPLE } from "${baseDir}/modules/plot_report/main"
//include { PLOT_REPORT_MULTI_SAMPLE } from "${baseDir}/modules/plot_report/plot_report_multi"

workflow TUMOUREVO {

  take:
  input_samplesheet
  fasta
  
  main:

  input = input_samplesheet.map{ meta, vcf, tbi, bam, bai, cna_segs, cna_extra -> 
        meta = meta + [id: "${meta.dataset}_${meta.patient}_${meta.tumour_sample}"]
        [ meta, vcf, tbi, bam, bai, cna_segs, cna_extra ]
      }


  input_vcf = input.map{ meta, vcf, tbi, bam, bai, cna_segs, cna_extra  -> 
        [ meta, vcf, tbi ]
  }

  input_cna = input.map{ meta, vcf, tbi, bam, bai, cna_segs, cna_extra  -> 
        [ meta, [cna_segs, cna_extra] ]
  }

  input_bam = input.map{ meta, vcf, tbi, bam, bai, cna_segs, cna_extra  -> 
        [ meta, bam, bai ]
  }


  ch_extra_files = []
  VCF_ANNOTATE_ENSEMBLVEP(input_vcf, 
                          fasta,
                          params.genome,
                          params.species,
                          params.vep_cache_version,
                          params.vep_dir_cache,
                          ch_extra_files)
  FORMATTER_VCF(VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi, "vcf") // return check multi-sample
  FORMATTER_CNA(input_cna, "cna")
  
  // lifter=false
  // if (params.mode == 'multisample' && lifter){ // + check vcf multi-sample
  // //  tumor_bam = Channel.fromPath(params.samples).
  // //    splitCsv(header: true).
  // //     map{row ->
  // //     tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.tumour_bam), file(row.tumour_bai))}

  //   LIFTER(FORMATTER_VCF.out, tumor_bam, fasta)
  //   annotation = DRIVER_ANNOTATION(LIFTER.out, cancer_type)

  // } else {
  //  annotation = DRIVER_ANNOTATION(FORMATTER_VCF.out, cancer_type)
  // }

  // QC(FORMATTER_CNA.out, annotation)
  // SUBCLONAL_DECONVOLUTION(QC.out.rds_join)
  // SIGNATURE_DECONVOLUTION(QC.out.rds_join)
}
