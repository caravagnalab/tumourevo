//
// MUTATIONAL SIGNATURES DECONVOLUTION WORKFLOW
//

include { FORMATTER as FORMATTER_RDS_SIGPROFILER } from "../../../subworkflows/local/formatter/main"
include { FORMATTER as FORMATTER_RDS_SPARSESIGNATURES } from "../../../subworkflows/local/formatter/main"
include { SPARSE_SIGNATURES } from "../../../modules/local/SparseSignatures/main"
include { DOWNLOAD_GENOME_SIGPROFILER } from "../../../modules/local/SigProfiler/download/main"
include { SIGPROFILER } from "../../../modules/local/SigProfiler/SigProfiler/main"


workflow SIGNATURE_DECONVOLUTION {
    take: 
    rds_join // tuple val(meta), path("*.rds"), val(tumor_samples)

    main:
    plot_pdf = null
    plot_rds = null
    signatures_nmfOut = null 
    bestConf = null
    sign_cv = null
    genome_path = null
    Sigprofiler_out = null


    if (params.tools && params.tools.split(',').contains('sparsesignatures')) {
        out_sparse = FORMATTER_RDS_SPARSESIGNATURES(rds_join, "rds")
        SPARSE_SIGNATURES(out_sparse.groupTuple(by: 0)) // run SparseSignatures
        
        plot_pdf = SPARSE_SIGNATURES.out.signatures_plot_pdf
        plot_rds = SPARSE_SIGNATURES.out.signatures_plot_rds
        signatures_nmfOut = SPARSE_SIGNATURES.out.signatures_nmfOut_rds
        bestConf = SPARSE_SIGNATURES.out.signatures_bestConf_rds
        sign_cv = SPARSE_SIGNATURES.out.signatures_cv_rds
    } 
  

    if (params.tools && params.tools.split(',').contains('sigprofiler')) {
        // Check if we should download SigProfiler genome
        if (params.download_sigprofiler_genome) {    
            genome_path = DOWNLOAD_GENOME_SIGPROFILER(params.download_genome_sigprofiler_reference_genome).genome_sigprofiler           
        } else {
        // Use the installed genome path from params
        genome_path = params.genome_installed_path
        }
        out_sigprof = FORMATTER_RDS_SIGPROFILER(rds_join, "rds")
        //Sigprofiler_out =  SIGPROFILER(out_sigprof.groupTuple(by: 0), genome_path) // run Sigprofiler 
        Sigprofiler_out =  SIGPROFILER(out_sigprof, genome_path)
    }

    emit:
    plot_pdf
    plot_rds
    signatures_nmfOut
    bestConf
    sign_cv
    Sigprofiler_out
    
}
